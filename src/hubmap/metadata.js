/**
 * @typedef {object} HubmapMetadata
 * @property {string} uuid
 * @property {string} organ
 * @property {string} organ_source
 * @property {string} assay_type
 * @property {string} dataset_iri
 * @property {string} donor_sex
 * @property {string} donor_race
 * @property {string} donor_age
 * @property {string} donor_bmi
 * @property {string} block_sample_iri
 * @property {string} section_sample_iri
 * @property {string} donor_iri
 */

const HUBMAP_ENTITY_ENDPOINT =
  'https://entity.api.hubmapconsortium.org/entities/';
const HUBMAP_PORTAL_ENDPOINT =
  'https://portal.hubmapconsortium.org/browse/dataset/';

const ORGAN_MAPPING = {
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

function getBody(ids) {
  return {
    version: true,
    from: 0,
    size: 10000,
    query: {
      terms: {
        'hubmap_id.keyword': ids,
      },
    },
    _source: {
      includes: [
        'uuid',
        'hubmap_id',
        'origin_samples.organ',
        'data_types',
        'mapped_consortium',
        'group_name',
        'group_uuid',
        'donor.mapped_metadata.race',
        'donor.mapped_metadata.sex',
        'donor.mapped_metadata.age_value',
        'donor.mapped_metadata.body_mass_index_value',
        'ancestors',
        'donor.uuid',
      ],
    },
  };
}

function checkResponse(response) {
  if (!response.ok) {
    const { status, statusText } = response;
    const message = `Failed to fetch metadata: ${status}:${statusText}`;
    throw new Error(message);
  }
}

function getSampleIri(ancestors, type) {
  for (const ancestor of ancestors) {
    if (ancestor['entity_type'].toLowerCase() == 'sample') {
      if (ancestor['sample_category'].toLowerCase() == type) {
        return `${HUBMAP_ENTITY_ENDPOINT}${ancestor['uuid']}`;
      }
    }
  }
  return '';
}

function toLookup(result) {
  /** @type {Map<string, HubmapMetadata>} */
  const lookup = new Map();
  for (const hit of result.hits.hits) {
    const {
      _source: {
        hubmap_id,
        uuid,
        origin_samples: [{ organ }],
        data_types: [assay_type],
        mapped_consortium,
        group_name,
        group_uuid,
        donor: {
          mapped_metadata: {
            age_value: [donor_age],
            race: [donor_race],
            sex: [donor_sex],
            body_mass_index_value: [donor_bmi],
          },
        },
        ancestors,
        donor: { uuid: donor_uuid },
      },
    } = hit;

    const mapped_organ = ORGAN_MAPPING[organ.toUpperCase()];

    lookup.set(hubmap_id, {
      organ: mapped_organ,
      organ_source: organ,
      uuid,
      assay_type,
      dataset_id: `${HUBMAP_ENTITY_ENDPOINT}${uuid}`,
      dataset_link: `${HUBMAP_PORTAL_ENDPOINT}${uuid}`,
      dataset_technology: 'OTHER',
      consortium_name: mapped_consortium,
      provider_name: group_name,
      provider_uuid: group_uuid,
      donor_id: `${HUBMAP_ENTITY_ENDPOINT}${donor_uuid}`,
      donor_age: donor_age ?? '',
      donor_sex: donor_sex ?? '',
      donor_bmi: donor_bmi ?? '',
      donor_race: donor_race ?? '',
      organ_id: `http://purl.obolibrary.org/obo/UBERON_${
        mapped_organ.split(':')[1]
      }`,
      block_id: getSampleIri(ancestors, 'block'),
      section_id: getSampleIri(ancestors, 'section'),
    });
  }

  return lookup;
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

export async function getMetadataLookup(ids, url, token) {
  let resp = await fetch(url, {
    method: 'POST',
    headers: getHeaders(token),
    body: JSON.stringify(getBody(ids)),
  });

  if (resp.status === 303) {
    resp = await handle303Response(resp);
  }

  checkResponse(resp);

  const result = await resp.json();
  return toLookup(result);
}
