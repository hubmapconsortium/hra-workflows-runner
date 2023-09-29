import { execFile as callbackExecFile } from 'node:child_process';
import { join } from 'node:path';
import { promisify } from 'node:util';
import { Cache } from '../util/cache.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { FORCE } from '../util/constants.js';
import { checkFetchResponse, downloadFile } from '../util/fs.js';
import { IDownloader } from '../util/handler.js';
import { getCacheDir, getSrcFilePath } from '../util/paths.js';
import { CollectionMetadata, parseMetadataFromId } from './metadata.js';
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
    /** @type {Cache<string, Promise<CollectionMetadata>>} */
    this.collectionCache = new Cache();
    /** @type {Cache<string, Promise<void>>} */
    this.assetCache = new Cache();
    /** @type {string} */
    this.extractScriptFile = 'extract_dataset.py';
    /** @type {string} */
    this.extractScriptFilePath = getSrcFilePath(
      config,
      'cellxgene',
      this.extractScriptFile
    );
  }

  async prepareDownload(datasets) {
    await concurrentMap(datasets, async (dataset) => {
      const metadata = parseMetadataFromId(dataset.id);
      const { dataset: datasetId, collection: collectionId, tissue } = metadata;
      if (collectionId) {
        const collection = await this.downloadCollection(collectionId);
        const { id: tissueId } = collection.findTissueId(datasetId, tissue);
        const { id: asset } = collection.findH5adAsset(datasetId) ?? {};
        Object.assign(dataset, metadata, { tissueId, asset });
      }
    });

    const tissueIds = datasets.map((dataset) => dataset.tissueId);
    const uniqueTissueIds = new Set(tissueIds);
    uniqueTissueIds.delete(undefined);

    const organLookup = await getOrganLookup(Array.from(uniqueTissueIds));
    for (const dataset of datasets) {
      dataset.organ = organLookup.get(dataset.tissueId);
    }

    return datasets.filter((dataset) => dataset.organ);
  }

  async download(dataset) {
    if (!dataset.asset) {
      throw new Error(`Cannot find h5ad data file`);
    }

    const assetDataFilePath = join(
      getCacheDir(this.config),
      `cellxgene-${dataset.asset}.h5ad`
    );

    await this.downloadAsset(dataset.dataset, dataset.asset, assetDataFilePath);
    await execFile('python3', [
      this.extractScriptFilePath,
      assetDataFilePath,
      '--donor',
      dataset.donor,
      '--tissue',
      dataset.tissue,
      '--sample',
      dataset.sample,
      '--output',
      dataset.dataFilePath,
    ]);
  }

  async downloadCollection(collection) {
    const url = new URL(`${COLLECTIONS_PATH}${collection}`, this.endpoint);
    return this.collectionCache.setDefaultFn(collection, () =>
      CollectionMetadata.download(url)
    );
  }

  async downloadAsset(dataset, asset, outputFile) {
    const pathname = `${ASSETS_PATH}${dataset}/asset/${asset}`;
    const url = new URL(pathname, this.endpoint);
    return this.assetCache.setDefaultFn(asset, async () => {
      const resp = await fetch(url, { method: 'POST' });
      checkFetchResponse(resp, 'CellXGene asset download failed');

      const { presigned_url } = await resp.json();
      await downloadFile(outputFile, presigned_url, {
        overwrite: this.config.get(FORCE, false),
      });
    });
  }
}
