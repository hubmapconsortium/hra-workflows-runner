import { join } from 'node:path';
import { Cache } from '../util/cache.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { checkFetchResponse, downloadFile } from '../util/fs.js';
import { IDownloader } from '../util/handler.js';
import { getCacheDir } from '../util/paths.js';
import { CollectionMetadata, parseMetadataFromId } from './metadata.js';

const CELLXGENE_API_ENDPOINT = 'CELLXGENE_API_ENDPOINT';
const DEFAULT_CELLXGENE_API_ENDPOINT = 'https://api.cellxgene.cziscience.com';
const COLLECTIONS_PATH = '/dp/v1/collections/';
const ASSETS_PATH = '/dp/v1/datasets/';

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
  }

  async prepareDownload(datasets) {
    await concurrentMap(datasets, async (dataset) => {
      const metadata = parseMetadataFromId(dataset.id);
      const { dataset: datasetId, collection: collectionId } = metadata;
      if (collectionId) {
        const collection = await this.downloadCollection(collectionId);
        const { id: asset } = collection.findH5adAsset(datasetId) ?? {};
        Object.assign(dataset, metadata, { asset });
      }
    });
  }

  async download(dataset) {
    if (!dataset.asset) {
      throw new Error(`Cannot find h5ad data file`);
    }

    const dataFilePath = join(
      getCacheDir(this.config),
      `cellxgene-${dataset.asset}.h5ad`
    );

    await this.downloadAsset(dataset.dataset, dataset.asset, dataFilePath);
    // TODO Use python script to split out from the h5ad file (like in gtex/download.js)
    // Should pass file, dataset.donor, dataset.tissue, and dataset.sample to the script
    // Script should use sample to split if a sample column exists otherwise it should
    // use donor combined with tissue to do the split
    // NOTE: Make sure donor, tissue, and sample is properly quoted when passed to the script
    // I.e. `... --tissue "tissue with spaces"` rather than `... --tissue tissue with spaces`
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
      await downloadFile(outputFile, presigned_url);
    });
  }
}
