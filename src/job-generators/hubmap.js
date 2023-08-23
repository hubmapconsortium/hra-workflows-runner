import { writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { env } from 'node:process';

const SEARCH_URL = 'https://search.api.hubmapconsortium.org/v3/portal/search';
const HUBMAP_TOKEN = env['HUBMAP_TOKEN'];

const RUI_ORGAN_MAPPING = {
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

async function organLookup(ids, token) {
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
      includes: ['hubmap_id', 'origin_samples.organ'],
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
  for (const item of result.hits.hits) {
    const {
      _source: {
        hubmap_id,
        origin_samples: [{ organ }],
      },
    } = item;
    mapping[hubmap_id] = RUI_ORGAN_MAPPING[organ];
  }

  return mapping;
}

function getJobYaml(organ) {
  return `organ: ${organ}
matrix:
  class: File
  path: data.h5ad

preprocessing:
  geneColumn: hugo_symbol

algorithms:
  - azimuth: {}
    extract:
      annotationMethod: azimuth
      cellLabelColumn: azimuth_label
      geneLabelColumn: gene

  - celltypist: {}
    extract:
      annotationMethod: celltypist
      cellLabelColumn: predicted_labels
      geneLabelColumn: gene

#  - popv:
#      referenceData:
#        class: File
#        path: TS_Lung_filtered.h5ad
#    extract:
#      annotationMethod: popv
#      cellLabelColumn: popv_prediction
#      geneLabelColumn: gene
`;
}

async function generateJob(dir, organ) {
  if (!organ) {
    return;
  }

  try {
    const path = join(dir.path, 'job.yaml');
    await writeFile(path, getJobYaml(organ));
  } catch (error) {
    console.error(error);
  }
}

export function supports(dir) {
  return /^hbm/i.test(dir.name);
}

export async function generate(dirs) {
  const ids = dirs.map((dir) => dir.name);
  const organMapping = await organLookup(ids, HUBMAP_TOKEN);
  for (const dir of dirs) {
    await generateJob(dir, organMapping[dir.name]);
  }
}
