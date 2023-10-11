import { execFile as callbackExecFile } from 'node:child_process';
import { writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { promisify } from 'node:util';
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
    /** @type {Cache<string, Promise<string>>} */
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
    const maxConcurrency = this.config.get(
      MAX_CONCURRENCY,
      DEFAULT_MAX_CONCURRENCY
    );
    const attachMetadata = async (dataset) => {
      const metadata = parseMetadataFromId(dataset.id);
      const { assets, tissueIdLookup } = await this.downloadCollection(
        metadata.collection
      );
      const tissue = metadata.tissue.toLowerCase();
      const tissueId = tissueIdLookup.get(tissue);
      Object.assign(dataset, metadata, { assets, tissueId });
    };

    await concurrentMap(datasets, attachMetadata, { maxConcurrency });
    datasets = datasets.filter(({ tissueId }) => !!tissueId);

    const tissueIds = datasets
      .map(({ tissueId }) => tissueId)
      .filter((id) => id.startsWith('UBERON:'))
      .map((id) => id.slice(0, 14));
    const uniqueTissueIds = Array.from(new Set(tissueIds));
    const organLookup = await getOrganLookup(uniqueTissueIds);
    for (const dataset of datasets) {
      dataset.organ = organLookup.get(dataset.tissueId);
    }

    datasets = datasets.filter(({ organ }) => !!organ);
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
    const errors = await this.extractDatasets(
      dataset.collection,
      dataset.assets,
      assetsFiles
    );

    this.checkDatasetExtractionError(dataset, errors);
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

  async extractDatasets(collection, assets, assetFiles) {
    return this.extractCache.setDefaultFn(collection, async () => {
      const datasets = this.datasetsByCollection.get(collection);
      const content = this.serializeExtractInfo(datasets, assets, assetFiles);
      const extractInfoFilePath = join(
        getCacheDir(this.config),
        `cellxgene-extract-info-${collection}.json`
      );

      await writeFile(extractInfoFilePath, content);
      const { stderr } = await execFile('python3', [
        this.extractScriptFilePath,
        extractInfoFilePath,
      ]);

      return stderr;
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
      throw new Error(msg);
    }
  }
}
