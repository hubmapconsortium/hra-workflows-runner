import { basename } from 'node:path';
import { checkFetchResponse } from '../util/fs.js';

export function parseMetadataFromId(id) {
  const { pathname, hash } = new URL(id);
  const collection = basename(pathname);
  const decodedHash = decodeURIComponent(hash);
  const [dataset, donor, tissue, sample = ''] = decodedHash.slice(1).split('$');
  return { collection, dataset, donor, tissue, sample };
}

export class CollectionMetadata {
  constructor(raw) {
    this.raw = raw;
  }

  findDataset(id) {
    const isPrimary = ({ is_primary_data: value }) => /primary/i.test(value);
    const isNonDiseased = ({ disease: [{ label }] }) => /normal/i.test(label);
    const isHuman = ({ organism: [{ label }] }) => /homo sapiens/i.test(label);
    return this.raw.datasets.find(
      (dataset) =>
        dataset.id === id &&
        isPrimary(dataset) &&
        isNonDiseased(dataset) &&
        isHuman(dataset)
    );
  }

  findH5adAsset(id) {
    const isH5ad = ({ filetype }) => /h5ad/i.test(filetype);
    const { dataset_assets: assets = [] } = this.findDataset(id) ?? [];
    return assets.find(isH5ad);
  }

  static async download(url) {
    const resp = await fetch(url, { method: 'GET' });
    checkFetchResponse(resp, 'CellXGene collection download failed');
    return new CollectionMetadata(await resp.json());
  }
}
