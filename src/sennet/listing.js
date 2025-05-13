import { IListing } from '../util/handler.js';
import { XConsortiaListing } from '../xconsortia/listing.js';

const SENNET_TOKEN = 'SENNET_TOKEN';
const SENNET_SEARCH_URL = 'SENNET_SEARCH_URL';
const DEFAULT_SENNET_SEARCH_URL = 'https://search.api.sennetconsortium.org/entities/search';

/** @implements {IListing} */
export class Listing extends XConsortiaListing {
  constructor(config) {
    super(config, config.get(SENNET_TOKEN), config.get(SENNET_SEARCH_URL, DEFAULT_SENNET_SEARCH_URL), 'sennet_id');
  }

  getBody() {
    return {
      version: true,
      from: 0,
      size: 10000,
      query: {
        bool: {
          must: [
            {
              term: {
                'entity_type.keyword': 'Dataset',
              },
            },
            {
              term: {
                'files.rel_path.keyword': 'expr.h5ad',
              },
            },
            {
              term: {
                'sources.source_type.keyword': 'Human',
              },
            },
          ],
        },
      },
      _source: {
        includes: [this.idKeyword],
      },
    };
  }
}
