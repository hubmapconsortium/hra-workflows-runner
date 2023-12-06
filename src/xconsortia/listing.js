import { Config } from '../util/config.js';
import { checkFetchResponse } from '../util/fs.js';
import { IListing } from '../util/handler.js';
import { getHeaders } from './metadata.js';

/** @implements {IListing} */
export class XConsortiaListing {
  constructor(config, token, searchUrl, idKeyword) {
    /** @type {Config} */
    this.config = config;
    /** @type {string} */
    this.token = token;
    /** @type {string} */
    this.searchUrl = searchUrl;
    /** @type {string} */
    this.idKeyword = idKeyword;
  }

  async getDatasets() {
    const resp = await fetch(this.searchUrl, {
      method: 'POST',
      headers: getHeaders(this.token),
      body: JSON.stringify(this.getBody()),
    });
    checkFetchResponse(resp, `${this.idKeyword.split('_')[0]}: Failed to fetch list of collections`);

    const {
      hits: { hits },
    } = await resp.json();
    return hits.map(({ _source }) => _source[this.idKeyword]);
  }

  getBody() {
    throw new Error('getBody must be overriden');
  }
}
