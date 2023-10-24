import { Config } from '../util/config.js';
import { IJobGenerator } from '../util/handler.js';

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

/** @implements {IJobGenerator} */
export class JobGenerator {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
  }

  async prepareJobs(datasets) {}

  createJob(dataset) {
    const organ = dataset.organ.toUpperCase();
    const uberon = ORGAN_MAPPING[organ];
    if (uberon === undefined) {
      const msg = `Unknown organ code '${dataset.organ}'`;
      throw new Error(msg);
    }

    return {
      organ: uberon,
      geneColumn: 'hugo_symbol',
      azimuth: {},
      celltypist: {},
      popv: {},
    };
  }
}
