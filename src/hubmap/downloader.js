import { OrganMetadataCollection } from '../organ/metadata.js';
import { IDownloader } from '../util/handler.js';
import { XConsortiaDownloader } from '../xconsortia/downloader.js';
import { getMetadata } from '../xconsortia/metadata.js';
import { ID_KEYWORD, METADATA_FIELDS, metadataToLookup } from './metadata.js';

const HUBMAP_TOKEN = 'HUBMAP_TOKEN';
const HUBMAP_SEARCH_URL = 'HUBMAP_SEARCH_URL';
const HUBMAP_ASSETS_URL = 'HUBMAP_ASSETS_URL';
const DEFAULT_HUBMAP_SEARCH_URL = 'https://search.api.hubmapconsortium.org/v3/portal/search';
const DEFAULT_HUBMAP_ASSETS_URL = 'https://assets.hubmapconsortium.org/';

/** @implements {IDownloader} */
export class Downloader extends XConsortiaDownloader {
  constructor(config) {
    super(
      config,
      config.get(HUBMAP_TOKEN),
      config.get(HUBMAP_SEARCH_URL, DEFAULT_HUBMAP_SEARCH_URL),
      config.get(HUBMAP_ASSETS_URL, DEFAULT_HUBMAP_ASSETS_URL)
    );
  }

  async getMetadataLookup(ids) {
    const organMetadata = await OrganMetadataCollection.load(this.config);
    const metadata = await getMetadata(ids, this.searchUrl, this.token, ID_KEYWORD, METADATA_FIELDS);
    return metadataToLookup(metadata, organMetadata);
  }
}
