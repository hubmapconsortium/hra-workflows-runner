import { Config } from '../util/config.js';
import { UnknownOrganError } from '../util/errors.js';
import { IJobGenerator } from '../util/handler.js';

const AUTHOR_CELL_TYPE_COLUMN = 'authors_cell_type_-_ontology_labels';

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
      geneColumn: 'feature_name',
      azimuth: {
        queryLayersKey: 'counts',
      },
      celltypist: {
        queryLayersKey: 'counts',
      },
      popv: {
        queryLayersKey: 'counts',
      },
      'pan-human-azimuth': {
        queryLayersKey: 'counts',
      },
      frmatch: {
        queryLayersKey: 'counts',
      },
      author: {
        labelColumn: AUTHOR_CELL_TYPE_COLUMN,
        organLevel: 'anatomogram',
        queryLayersKey: 'counts',
      },
      cellSource: dataset.dataset_id,
    };
  }
}
