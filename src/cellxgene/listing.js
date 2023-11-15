import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { DEFAULT_MAX_CONCURRENCY, MAX_CONCURRENCY } from '../util/constants.js';
import { checkFetchResponse } from '../util/fs.js';
import { IListing } from '../util/handler.js';
import { downloadCollectionMetadata } from './metadata.js';

const CELLXGENE_API_ENDPOINT = 'CELLXGENE_API_ENDPOINT';
const DEFAULT_CELLXGENE_API_ENDPOINT = 'https://api.cellxgene.cziscience.com';
const COLLECTIONS_PATH = '/dp/v1/collections/';

/** @implements {IListing} */
export class Listing {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.endpoint = config.get(
      CELLXGENE_API_ENDPOINT,
      DEFAULT_CELLXGENE_API_ENDPOINT
    );
  }

  async getDatasets() {
    const maxConcurrency = this.config.get(
      MAX_CONCURRENCY,
      DEFAULT_MAX_CONCURRENCY
    );

    const collectionIds = await this.getCollections();
    const urls = collectionIds.map(
      (id) => new URL(`${COLLECTIONS_PATH}${id}`, this.endpoint)
    );
    const metadata = await concurrentMap(urls, downloadCollectionMetadata, {
      maxConcurrency,
    });

    return metadata.flatMap(({ id, donorTissuePairs }) =>
      donorTissuePairs.map(({ donor_id, tissue }) => {
        const encodedDonor = encodeURIComponent(donor_id);
        const encodedTissue = encodeURIComponent(tissue);
        return `${this.endpoint}${COLLECTIONS_PATH}${id}#${encodedDonor}$${encodedTissue}`;
      })
    );
  }

  async getCollections() {
    const url = new URL(COLLECTIONS_PATH, this.endpoint);
    const resp = await fetch(url, { method: 'GET' });
    checkFetchResponse(resp, 'CellXGene: Failed to fetch list of collections');

    const { collections } = await resp.json();
    return collections.map(({ id }) => id);
  }
}
