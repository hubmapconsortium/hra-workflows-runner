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

function toLookup(result) {
  /** @type {Map<string, { uuid: string; organ: string; assay_type: string; dataset_iri: string, donor_sex: string, donor_race: string, donor_age: string}>} */
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
      },
    } = hit;
    const HUBMAP_ENTITY_ENDPOINT =
      'https://entity.api.hubmapconsortium.org/entities/';
    const dataset_iri = HUBMAP_ENTITY_ENDPOINT.concat(uuid);
    const donor_age = age_value + ' '.concat(age_unit);
    lookup.set(hubmap_id, {
      uuid,
      organ,
      assay_type,
      dataset_iri,
      donor_sex,
      donor_race,
      donor_age,
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
