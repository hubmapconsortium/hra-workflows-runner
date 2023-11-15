import { Config } from '../util/config.js';
import { checkFetchResponse } from '../util/fs.js';
import { IListing } from '../util/handler.js';
import { getHeaders } from './metadata.js';

const HUBMAP_TOKEN = 'HUBMAP_TOKEN';
const HUBMAP_SEARCH_URL = 'HUBMAP_SEARCH_URL';
const DEFAULT_HUBMAP_SEARCH_URL =
  'https://search.api.hubmapconsortium.org/v3/portal/search';

/** @implements {IListing} */
export class Listing {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.token = config.get(HUBMAP_TOKEN);
    /** @type {string} */
    this.searchUrl = config.get(HUBMAP_SEARCH_URL, DEFAULT_HUBMAP_SEARCH_URL);
  }

  async getDatasets() {
    const resp = await fetch(this.searchUrl, {
      method: 'POST',
      headers: getHeaders(this.token),
      body: JSON.stringify(this.getBody()),
    });
    checkFetchResponse(resp, 'HuBMAP: Failed to fetch list of collections');

    const {
      hits: { hits },
    } = await resp.json();
    return hits.map(({ _source: { hubmap_id } }) => hubmap_id);
  }

  getBody() {
    return {
      version: true,
      from: 0,
      size: 10000,
      query: {
        term: {
          'files.rel_path.keyword': 'expr.h5ad',
        },
      },
      _source: {
        includes: ['hubmap_id'],
      },
    };
  }
}
