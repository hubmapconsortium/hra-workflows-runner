import { Config } from '../util/config.js';
import { IJobGenerator } from '../util/handler.js';

const ORGAN_MAPPING = {
  bladder: 'UBERON:0001255',
  blood: 'UBERON:0000178',
  bone_marrow: 'UBERON:0002371',
  eye: 'UBERON:0000970',
  heart: 'UBERON:0000948',
  large_intestine: 'UBERON:0000059',
  liver: 'UBERON:0002107',
  lung: 'UBERON:0002048',
  lymph_node: 'UBERON:0000029',
  mammary: 'UBERON:0001911',
  pancreas: 'UBERON:0001264',
  prostate: 'UBERON:0002367',
  skin: 'UBERON:0002097',
  small_intestine: 'UBERON:0002108',
  spleen: 'UBERON:0002106',
  thymus: 'UBERON:0002370',
  trachea: 'UBERON:0003126',
  uterus: 'UBERON:0000995',
  vasculature: 'UBERON:0004537',
};

/** @implements {IJobGenerator} */
export class JobGenerator {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
  }

  async prepareJobs(datasets) {}

  createJob(dataset) {
    const organ = dataset.organ.toLowerCase();
    const uberon = ORGAN_MAPPING[organ];
    if (uberon === undefined) {
      const msg = `Unknown organ code '${dataset.organ}'`;
      throw new Error(msg);
    }

    return {
      organ: uberon,
      geneColumn: 'gene_name',
      azimuth: {},
      celltypist: {},
      popv: {},
    };
  }
}
