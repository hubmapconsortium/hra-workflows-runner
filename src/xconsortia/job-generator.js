import { Config } from '../util/config.js';
import { UnknownOrganError } from '../util/errors.js';
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
      throw new UnknownOrganError(dataset.organ_source);
    }

    return {
      organ: dataset.organ,
      geneColumn: 'hugo_symbol',
      azimuth: {
        queryLayersKey: 'raw',
      },
      celltypist: {
        queryLayersKey: 'raw',
      },
      popv: {
        queryLayersKey: 'raw',
      },
      'pan-human-azimuth': {
        queryLayersKey: 'raw',
      },
      'frmatch': {
        queryLayersKey: 'raw', 
      },
      cellSource: dataset.dataset_id,
    };
  }
}
