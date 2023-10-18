import { spawn } from 'node:child_process';
import { writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { Dataset } from '../dataset/dataset.js';
import { Cache } from '../util/cache.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import {
  DEFAULT_MAX_CONCURRENCY,
  FORCE,
  MAX_CONCURRENCY,
} from '../util/constants.js';
import { checkFetchResponse, downloadFile } from '../util/fs.js';
import { IDownloader } from '../util/handler.js';
import { groupBy } from '../util/iter.js';
import { getCacheDir, getSrcFilePath } from '../util/paths.js';
import { downloadCollectionMetadata, parseMetadataFromId } from './metadata.js';
import { getOrganLookup } from './organ-lookup.js';

const CELLXGENE_API_ENDPOINT = 'CELLXGENE_API_ENDPOINT';
const DEFAULT_CELLXGENE_API_ENDPOINT = 'https://api.cellxgene.cziscience.com';
const COLLECTIONS_PATH = '/dp/v1/collections/';
const ASSETS_PATH = '/dp/v1/datasets/';
const UBERON_ID_REGEX = /^UBERON:\d{7}$/;

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

    const error = dataset.scratch.extractError;
    if (error) {
      throw error;
    }
  }

  async downloadCollection(collection) {
    const url = new URL(`${COLLECTIONS_PATH}${collection}`, this.endpoint);
    return this.collectionCache.setDefaultFn(collection, () =>
      downloadCollectionMetadata(url)
    );
  }

  async downloadAssets(/** @type {any[]} */ assets) {
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

        await downloadFile(outputFile, presigned_url, {
          overwrite: this.config.get(FORCE, false),
        });

        return outputFile;
      });

    const downloads = assets.map(downloadAsset);
    return Promise.all(downloads);

    // const pathname = `${ASSETS_PATH}${dataset}/asset/${asset}`;
    // const url = new URL(pathname, this.endpoint);
    // return this.assetCache.setDefaultFn(asset, async () => {
    //   const resp = await fetch(url, { method: 'POST' });
    //   checkFetchResponse(resp, 'CellXGene asset download failed');

    //   const { presigned_url } = await resp.json();
    //   await downloadFile(outputFile, presigned_url, {
    //     overwrite: this.config.get(FORCE, false),
    //   });
    // });
  }

  async attachMetadata(datasets) {
    const maxConcurrency = this.config.get(
      MAX_CONCURRENCY,
      DEFAULT_MAX_CONCURRENCY
    );
    const attach = async (dataset) => {
      const metadata = parseMetadataFromId(dataset.id);
      const { assets, tissueIdLookup } = await this.downloadCollection(
        metadata.collection
      );
      const tissue = metadata.tissue.toLowerCase();
      const tissueId = tissueIdLookup.get(tissue);
      Object.assign(dataset, metadata, { assets, tissueId });
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
      try {
        console.log('CellXGene:Extract:Start', collection);
        await this.runExtractDatasetsScript(
          datasets,
          extractInfoFilePath,
          tempExtractDirPath
        );
      } finally {
        console.log('CellXGene:Extract:End', collection);
      }
    });
  }

  async runExtractDatasetsScript(
    datasets,
    extractInfoFilePath,
    tempExtractDirPath
  ) {
    return new Promise((resolve, reject) => {
      const process = spawn(
        'python3',
        [
          this.extractScriptFilePath,
          extractInfoFilePath,
          '--tmp-dir',
          tempExtractDirPath,
        ],
        { stdio: [null, 'inherit', 'pipe'] }
      );

      process.stderr.on('data', (data) => {
        for (const dataset of datasets) {
          this.checkDatasetExtractionError(dataset, data);
        }
      });

      process.on('error', reject);
      process.on('close', resolve);
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

  checkDatasetExtractionError(dataset, errors) {
    const { id } = dataset;
    const { length } = id;
    const start = errors.indexOf(id);
    if (start >= 0) {
      const end = errors.indexOf('\n', start);
      const msg = errors.slice(start + length + 2, end);
      dataset.scratch.extractError = new Error(msg);
    }
  }
}
