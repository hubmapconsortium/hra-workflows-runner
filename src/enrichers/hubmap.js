import { env } from 'node:process';

const SEARCH_URL = 'https://search.api.hubmapconsortium.org/v3/portal/search';
const HUBMAP_TOKEN = env['HUBMAP_TOKEN'];

/**
 * Looks up hubmap uuid for hubmap identifiers
 * @param {string[]} ids Identifiers to lookup
 * @param {?string} token Hubmap token
 * @returns {Promise<Object.<string, string>>} Mapping from hubmap ids to uuids
 */
async function uuidLookup(ids, token) {
  const headers = { 'Content-type': 'application/json' };
  const body = {
    version: true,
    from: 0,
    size: 10000,
    query: {
      terms: {
        'hubmap_id.keyword': ids,
      },
    },
    _source: {
      includes: ['uuid', 'hubmap_id'],
    },
  };

  if (token) {
    headers['Authorization'] = `Bearer ${token}`;
  }

  const resp = await fetch(SEARCH_URL, {
    method: 'POST',
    headers,
    body: JSON.stringify(body),
  });
  if (!resp.ok) {
    return {};
  }

  const result = await resp.json();
  const mapping = {};
  for (const { _source: { hubmap_id, uuid } } of result.hits.hits) {
    mapping[hubmap_id] = uuid;
  }

  return mapping;
}

/**
 * Enrich hubmap datasets with uuids
 * @param {object[]} datasets Datasets to enrich
 */
export default async function enrichHubmapDatasets(datasets) {
  const hubmapDatasets = datasets.filter(({ id }) => /^hbm/i.test(id));
  const ids = hubmapDatasets.map(({ id }) => id);
  const mapping = await uuidLookup(ids, HUBMAP_TOKEN);
  for (const dataset of hubmapDatasets) {
    dataset.uuid = mapping[dataset.id];
  }
}
