import { Config } from '../util/config.js';
import { IListing } from '../util/handler.js';

/** @implements {IListing} */
export class Listing {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
  }

  async getDatasets() {
    //
    return [];
  }
}
