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
      geneColumn: 'gene_symbol', // 'ensembl_id' or 'gene_symbol',
      azimuth: {
        queryLayersKey: 'raw_counts',
      },
      celltypist: {
        queryLayersKey: 'raw_counts',
      },
      popv: {
        queryLayersKey: 'raw_counts',
      },
      cellSource: dataset.dataset_id,
    };
  }
}
