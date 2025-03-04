import { IListing } from '../util/handler.js';
import { XConsortiaListing } from '../xconsortia/listing.js';

const HUBMAP_TOKEN = 'HUBMAP_TOKEN';
const HUBMAP_SEARCH_URL = 'HUBMAP_SEARCH_URL';
const DEFAULT_HUBMAP_SEARCH_URL = 'https://search.api.hubmapconsortium.org/v3/portal/search';

/** @implements {IListing} */
export class Listing extends XConsortiaListing {
  constructor(config) {
    super(config, config.get(HUBMAP_TOKEN), config.get(HUBMAP_SEARCH_URL, DEFAULT_HUBMAP_SEARCH_URL), 'hubmap_id');
  }

  getBody() {
    return {
      version: true,
      from: 0,
      size: 10000,
      query: {
        terms: {
          'files.rel_path.keyword': ['expr.h5ad'], // , 'cell_by_gene.h5ad'],
        },
      },
      _source: {
        includes: [this.idKeyword],
      },
    };
  }
}
