import { execFile as callbackExecFile } from 'node:child_process';
import { promisify } from 'node:util';
import { spawn } from 'node:child_process';
import { writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { Dataset } from '../dataset/dataset.js';
import { Cache } from '../util/cache.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import {
  DEFAULT_MAX_CONCURRENCY,
  DEFAULT_PYTHON_LOG_LEVEL,
  FORCE,
  MAX_CONCURRENCY,
  PYTHON_LOG_LEVEL,
} from '../util/constants.js';
import { checkFetchResponse, downloadFile, fileExists } from '../util/fs.js';
import { IDownloader } from '../util/handler.js';
import { groupBy } from '../util/iter.js';
import { logEvent } from '../util/logging.js';
import { getCacheDir, getSrcFilePath } from '../util/paths.js';
import { downloadCollectionMetadata, parseMetadataFromId } from './metadata.js';
import { getOrganLookup } from './organ-lookup.js';

const CELLXGENE_API_ENDPOINT = 'CELLXGENE_API_ENDPOINT';
const DEFAULT_CELLXGENE_API_ENDPOINT = 'https://api.cellxgene.cziscience.com';
const COLLECTIONS_PATH = '/dp/v1/collections/';
const ASSETS_PATH = '/dp/v1/datasets/';
const UBERON_ID_REGEX = /^UBERON:\d{7}$/;

const execFile = promisify(callbackExecFile);

/** @implements {IDownloader} */
export class Downloader {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.endpoint = config.get(
      CELLXGENE_API_ENDPOINT,
      DEFAULT_CELLXGENE_API_ENDPOINT
    );
    /** @type {Map<string, Dataset[]>} */
    this.datasetsByCollection = new Map();
    /** @type {Cache<string, Promise<import('./metadata.js').CollectionMetadata>>} */
    this.collectionCache = new Cache();
    /** @type {Cache<string, Promise<string>>} */
    this.assetCache = new Cache();
    /** @type {Cache<string, Promise<void>>} */
    this.extractCache = new Cache();
    /** @type {string} */
    this.extractScriptFile = 'extract_dataset_multi.py';
    /** @type {string} */
    this.extractScriptFilePath = getSrcFilePath(
      config,
      'cellxgene',
      this.extractScriptFile
    );
    this.extractMetdataScriptFile = 'extract_donor_metadata.py';
    /** @type {string} */
    this.extractMetdataScriptFilePath = getSrcFilePath(
      config,
      'cellxgene',
      this.extractMetdataScriptFile
    );
  }

  async prepareDownload(datasets) {
    datasets = await this.attachMetadata(datasets);
    datasets = await this.lookupOrgan(datasets);
    this.datasetsByCollection = groupBy(
      datasets,
      ({ collection }) => collection
    );

    return datasets;
  }

  async download(dataset) {
    if (dataset.assets.length === 0) {
      throw new Error(`Cannot find h5ad data file`);
    }

    const assetsFiles = await this.downloadAssets(dataset.assets);
    await this.extractDatasets(dataset.collection, dataset.assets, assetsFiles);
    if (!(await fileExists(dataset.dataFilePath))) {
      throw new Error('No data for donor-tissue combination');
    }

    const { stdout } = await execFile('python3', [
      this.extractMetdataScriptFilePath,
      dataset.dataFilePath,
    ]);

    // Parse sex line. Format: `sex: X\n`
    const sex_match = /sex:(.+)\n/i.exec(stdout);
    dataset.donor_sex = sex_match?.[1].trim() ?? '';

    // Parse age line. Format: `age: X\n`
    const age_match = /age:(.+)\n/i.exec(stdout);
    dataset.donor_age = age_match?.[1].trim() ?? '';

    // Parse age line. Format: `age: X\n`
    const ethnicity_match = /ethnicity:(.+)\n/i.exec(stdout);
    dataset.donor_race = ethnicity_match?.[1].trim() ?? '';
  }

  async downloadCollection(collection) {
    const url = new URL(`${COLLECTIONS_PATH}${collection}`, this.endpoint);
    return this.collectionCache.setDefaultFn(collection, () =>
      logEvent('CellXGene:DownloadCollection', collection, () =>
        downloadCollectionMetadata(url)
      )
    );
  }

  async downloadAssets(/** @type {any[]} */ assets) {
    const maxConcurrency = this.config.get(
      MAX_CONCURRENCY,
      DEFAULT_MAX_CONCURRENCY
    );
    const downloadAsset = ({ id, dataset }) =>
      this.assetCache.setDefaultFn(id, async () => {
        const pathname = `${ASSETS_PATH}${dataset}/asset/${id}`;
        const url = new URL(pathname, this.endpoint);
        const resp = await fetch(url, { method: 'POST' });
        checkFetchResponse(resp, 'CellXGene asset download failed');

        const { presigned_url } = await resp.json();
        const outputFile = join(
          getCacheDir(this.config),
          `cellxgene-${id}.h5ad`
        );

        await logEvent('CellXGene:DownloadAsset', id, () =>
          downloadFile(outputFile, presigned_url, {
            overwrite: this.config.get(FORCE, false),
          })
        );

        return outputFile;
      });

    return concurrentMap(assets, downloadAsset, {
      maxConcurrency,
    });
  }

  async attachMetadata(datasets) {
    const maxConcurrency = this.config.get(
      MAX_CONCURRENCY,
      DEFAULT_MAX_CONCURRENCY
    );
    const attach = async (dataset) => {
      const metadata = parseMetadataFromId(dataset.id);
      const { assets, tissueIdLookup, publicationDOI, donorTissuePairs } =
        await this.downloadCollection(metadata.collection);
      const tissue = metadata.tissue.toLowerCase();
      const tissueId = tissueIdLookup.get(tissue);
      Object.assign(dataset, metadata, {
        assets,
        tissueId,
        publicationDOI,
        donorTissuePairs,
      });
    };

    await concurrentMap(datasets, attach, { maxConcurrency });
    return datasets.filter(({ tissueId }) => !!tissueId);
  }

  async lookupOrgan(datasets) {
    const tissueIds = datasets
      .map(({ tissueId }) => tissueId)
      .filter((id) => UBERON_ID_REGEX.test(id));
    const uniqueTissueIds = Array.from(new Set(tissueIds));
    const organLookup = await getOrganLookup(uniqueTissueIds);

    for (const dataset of datasets) {
      dataset.organ = organLookup.get(dataset.tissueId) ?? '';
      if (dataset.organ === '') {
        const msg = `Cannot determine organ for tissue '${dataset.tissue}' (${dataset.tissueId})`;
        dataset.scratch.summary_ref.comments = msg;
      }
    }

    return datasets;
  }

  async extractDatasets(collection, assets, assetFiles) {
    return this.extractCache.setDefaultFn(collection, async () => {
      const datasets = this.datasetsByCollection.get(collection);
      const content = this.serializeExtractInfo(datasets, assets, assetFiles);
      const extractInfoFilePath = join(
        getCacheDir(this.config),
        `cellxgene-extract-info-${collection}.json`
      );
      const tempExtractDirPath = join(
        getCacheDir(this.config),
        `cellxgene-extract-${collection}`
      );

      await writeFile(extractInfoFilePath, content);
      await logEvent('CellXGene:Extract', collection, () =>
        this.runExtractDatasetsScript(extractInfoFilePath, tempExtractDirPath)
      );
    });
  }

  async runExtractDatasetsScript(extractInfoFilePath, tempExtractDirPath) {
    return new Promise((resolve, reject) => {
      const process = spawn(
        'python3',
        [
          this.extractScriptFilePath,
          extractInfoFilePath,
          '--tmp-dir',
          tempExtractDirPath,
          '--log-level',
          this.config.get(PYTHON_LOG_LEVEL, DEFAULT_PYTHON_LOG_LEVEL),
        ],
        { stdio: [null, 'inherit', 'inherit'] }
      );

      process.on('error', reject);
      process.on('exit', (code) => {
        if (code === 0) {
          resolve();
        } else {
          const msg = `Python exited with code ${code}`;
          reject(new Error(msg));
        }
      });
    });
  }

  serializeExtractInfo(datasets, assets, assetFiles) {
    const content = {
      assets,
      assetFiles,
      datasets: datasets.map((dataset) => ({
        id: dataset.id,
        donor: dataset.donor,
        tissue: dataset.tissue,
        outputFile: dataset.dataFilePath,
      })),
    };

    return JSON.stringify(content, undefined, 4);
  }
}
