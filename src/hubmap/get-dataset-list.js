import { getConfig } from '../util/common.js';

const HUBMAP_TOKEN =
  'AgqkKmq0JN50EpK02r1qklja26aPk7N54Jz0oYm243zN6njoW3UqCpWlNw2lldEYJ7mG7YVd5gXvmC5gjvYvtGd28k';
// const HUBMAP_SEARCH_URL = 'HUBMAP_SEARCH_URL';
const DEFAULT_HUBMAP_SEARCH_URL =
  'https://search.api.hubmapconsortium.org/v3/portal/search';

// https://portal.hubmapconsortium.org/browse/dataset/b6adb50fe8ab0d3d596a04b233d5c702.json

async function downloadAllDatasetListing() {
  const token = HUBMAP_TOKEN;
  const searchUrl = DEFAULT_HUBMAP_SEARCH_URL;

  fetch(searchUrl, {
    method: 'POST',
    headers: token
      ? { 'Content-type': 'application/json', Authorization: `Bearer ${token}` }
      : { 'Content-type': 'application/json' },
    body: JSON.stringify({
      version: true,
      from: 0,
      size: 10000,
      query: {
        term: {
          'files.rel_path.keyword': 'raw_expr.h5ad',
        },
      },
      _source: {
        includes: ['uuid', 'hubmap_id', 'data_types'],
      },
    }),
  })
    .catch((r) => ({ ok: false }))
    .then((r) => (r.ok ? r.json() : undefined))
    .then((r) =>
      r.hits.hits.map((hit) => ({
        ...hit._source,
        data_types: undefined,
        assay_type: hit._source.data_types.join(','),
      }))
    )
    .then((r) => console.log(JSON.stringify(r, null, 2), r.length));
}

async function main() {
  //   const config = getConfig();
  //   console.log(config);
  downloadAllDatasetListing();
}

main();
