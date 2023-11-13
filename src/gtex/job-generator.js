import { Config } from '../util/config.js';
import { IJobGenerator } from '../util/handler.js';

/** @implements {IJobGenerator} */
export class JobGenerator {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
  }

  async prepareJobs(datasets) {}

  createJob(dataset) {
    if (!dataset.organ) {
      const msg = `Unknown organ code '${dataset.organ_source}'`;
      throw new Error(msg);
    }

    return {
      organ: dataset.organ,
      geneColumn: 'gene_name',
      azimuth: {},
      celltypist: {},
      popv: {
        queryLayersKey: 'counts',
      },
    };
  }
}
