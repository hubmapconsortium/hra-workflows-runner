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
    return {
      organ: dataset.organ,
      geneColumn: 'feature_name',
      azimuth: {},
      celltypist: {},
      popv: false,
    };
  }
}
