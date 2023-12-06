import { IDownloader } from '../util/handler.js';
import { XConsortiaDownloader } from '../xconsortia/downloader.js';
import { getMetadata } from '../xconsortia/metadata.js';
import { ID_KEYWORD, METADATA_FIELDS, toLookup } from './metadata.js';

const SENNET_TOKEN = 'SENNET_TOKEN';
const SENNET_SEARCH_URL = 'SENNET_SEARCH_URL';
const SENNET_ASSETS_URL = 'SENNET_ASSETS_URL';
const DEFAULT_SENNET_SEARCH_URL = 'https://search.api.sennetconsortium.org/entities/search';
const DEFAULT_SENNET_ASSETS_URL = 'https://assets.api.sennetconsortium.org/';

/** @implements {IDownloader} */
export class Downloader extends XConsortiaDownloader {
  constructor(config) {
    super(
      config,
      config.get(SENNET_TOKEN),
      config.get(SENNET_SEARCH_URL, DEFAULT_SENNET_SEARCH_URL),
      config.get(SENNET_ASSETS_URL, DEFAULT_SENNET_ASSETS_URL)
    );
  }

  async getMetadataLookup(ids) {
    const metadata = await getMetadata(ids, this.searchUrl, this.token, ID_KEYWORD, METADATA_FIELDS);
    return toLookup(metadata);
  }
}
