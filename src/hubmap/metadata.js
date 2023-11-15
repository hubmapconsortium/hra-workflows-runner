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
 */

const HUBMAP_ENTITY_ENDPOINT =
  'https://entity.api.hubmapconsortium.org/entities/';

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
        'donor.mapped_metadata.race',
        'donor.mapped_metadata.sex',
        'donor.mapped_metadata.age_value',
        'donor.mapped_metadata.age_unit',
        'ancestors',
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
        donor: {
          mapped_metadata: {
            sex: [donor_sex],
            race: [donor_race],
            age_value: [age_value],
            age_unit: [age_unit],
          },
        },
        ancestors,
      },
    } = hit;

    lookup.set(hubmap_id, {
      uuid,
      organ_source: organ,
      organ: ORGAN_MAPPING[organ.toUpperCase()],
      assay_type,
      dataset_iri: `${HUBMAP_ENTITY_ENDPOINT}${uuid}`,
      donor_sex,
      donor_race,
      donor_age: `${age_value} ${age_unit}`,
      sample_block_iri: getSampleIri(ancestors, 'block'),
      sample_section_iri: getSampleIri(ancestors, 'section'),
    });
  }

  return lookup;
}

export async function getMetadataLookup(ids, url, token) {
  const resp = await fetch(url, {
    method: 'POST',
    headers: getHeaders(token),
    body: JSON.stringify(getBody(ids)),
  });
  checkResponse(resp);

  const result = await resp.json();
  return toLookup(result);
}
