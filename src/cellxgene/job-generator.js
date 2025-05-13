import { stat } from 'fs/promises';
import { Config } from '../util/config.js';
import { DATASET_DATA_MAX_SIZE, DEFAULT_DATASET_DATA_MAX_SIZE } from '../util/constants.js';
import { DataTooLargeError, UnknownOrganError } from '../util/errors.js';
import { IJobGenerator } from '../util/handler.js';

/** @implements {IJobGenerator} */
export class JobGenerator {
  constructor(config) {
    /** @type {Config} */
    this.config = config;
  }

  async prepareJobs(datasets) {}

  async createJob(dataset) {
    if (!dataset.organ) {
      throw new UnknownOrganError(`${dataset.tissue} (${dataset.tissueId})`);
    }

    const { size } = await stat(dataset.dataFilePath);
    const maxSize = this.config.get(DATASET_DATA_MAX_SIZE, DEFAULT_DATASET_DATA_MAX_SIZE);
    if (size > maxSize) {
      throw new DataTooLargeError(size);
    }

    return {
      organ: dataset.organ,
      geneColumn: 'feature_name',
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
      cellSource: dataset.dataset_id,
    };
  }
}
