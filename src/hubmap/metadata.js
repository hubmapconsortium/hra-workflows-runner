/**
 * @typedef {object} HubmapMetadata
 * @property {string} uuid
 * @property {string} organ
 * @property {string} assay_type
 * @property {string} dataset_iri
 * @property {string} donor_sex
 * @property {string} donor_race
 * @property {string} donor_age
 */

const HUBMAP_ENTITY_ENDPOINT =
  'https://entity.api.hubmapconsortium.org/entities/';

function getHeaders(token) {
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
      organ,
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
