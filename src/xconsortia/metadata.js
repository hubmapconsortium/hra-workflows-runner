import { checkFetchResponse } from '../util/fs.js';

export const ORGAN_MAPPING = {
  AO: 'UBERON:0000948',
  BL: 'UBERON:0001255',
  BD: 'UBERON:0001270',
  BM: 'UBERON:0001270',
  BR: 'UBERON:0000955',
  LB: 'UBERON:0001004',
  RB: 'UBERON:0001004',
  LE: 'UBERON:0004548',
  RE: 'UBERON:0004549',
  LF: 'UBERON:0001303',
  RF: 'UBERON:0001302',
  HT: 'UBERON:0000948',
  LK: 'UBERON:0004538',
  RK: 'UBERON:0004539',
  LI: 'UBERON:0000059',
  LV: 'UBERON:0002107',
  LL: 'UBERON:0001004',
  LN: 'FMA:24978',
  RL: 'UBERON:0001004',
  RN: 'FMA:24977',
  LY: 'UBERON:0002509',
  LO: 'FMA:7214',
  RO: 'FMA:7213',
  PA: 'UBERON:0001264',
  PL: 'UBERON:0001987',
  SI: 'UBERON:0002108',
  SK: 'UBERON:0002097',
  SP: 'UBERON:0002106',
  TH: 'UBERON:0002370',
  TR: 'UBERON:0001004',
  UR: 'UBERON:0001223',
  UT: 'UBERON:0000995',
};

export function getHeaders(token) {
  return {
    'Content-type': 'application/json',
    Authorization: `Bearer ${token}`,
  };
}

function getBody(ids, id_keyword, fields) {
  return {
    version: true,
    from: 0,
    size: 10000,
    query: {
      terms: {
        [`${id_keyword}.keyword`]: ids,
      },
    },
    _source: {
      includes: fields,
    },
  };
}

/**
 * Handles 303 responses from the search api.
 * A 303 response is returned when the resulting query is to large for the
 * search api. Instead it returns a temporary url from which to download the result.
 *
 * @param {Response} resp
 * @returns {Promise<Response>}
 */
async function handle303Response(resp) {
  const text = await resp.text();
  if (text.startsWith('https')) {
    return await fetch(text);
  }

  return resp;
}

export async function getMetadata(ids, url, token, id_keyword, fields) {
  let resp = await fetch(url, {
    method: 'POST',
    headers: getHeaders(token),
    body: JSON.stringify(getBody(ids, id_keyword, fields)),
  });
  if (resp.status === 303) {
    resp = await handle303Response(resp);
  }

  checkFetchResponse(resp, 'Failed to fetch metadata');
  const result = await resp.json();
  return result;
}

export function getSampleBlockId(ancestors, url_prefix) {
  for (const ancestor of ancestors) {
    if (ancestor['entity_type'].toLowerCase() == 'sample' && ancestor['sample_category'].toLowerCase() == 'block') {
      return {
        block_id: `${url_prefix}${ancestor['uuid']}`,
        rui_location: ancestor['rui_location'] ?? '',
      };
    }
  }

  return {};
}

export function getSampleSectionId(ancestors, url_prefix) {
  for (const ancestor of ancestors) {
    if (ancestor['entity_type'].toLowerCase() == 'sample' && ancestor['sample_category'].toLowerCase() == 'section') {
      return `${url_prefix}${ancestor['uuid']}`;
    }
  }

  return '';
}
