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
      includes: ['uuid', 'hubmap_id', 'origin_samples.organ', 'data_types'],
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
  /** @type {Map<string, { uuid: string; organ: string; assay_type: string; }>} */
  const lookup = new Map();
  for (const hit of result.hits.hits) {
    const {
      _source: {
        hubmap_id,
        uuid,
        origin_samples: [{ organ }],
        data_types: [ assay_type ],
      },
    } = hit;

    lookup.set(hubmap_id, { uuid, organ, assay_type });
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
