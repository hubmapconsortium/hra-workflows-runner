import { readFile, writeFile } from 'node:fs/promises';
import Papa from 'papaparse';
import { ALGORITHMS } from '../util/constants.js';

/**
 * Status enum for different dataset steps
 * @readonly
 * @enum {string}
 */
export const Status = {
  SUCCESS: 'success',
  FAILURE: 'failure',
  NOT_STARTED: 'not started',
  SKIPPED: 'skipped',
};

export class DatasetSummary {
  static fromRaw(data) {
    return Object.assign(new DatasetSummary(data.id), data);
  }

  /**
   * Creates a new dataset summary item
   *
   * @param {string} id Dataset id
   */
  constructor(id) {
    /** @type {string} */
    this.id = id;
    /** @type {Status} */
    this.downloaded = Status.NOT_STARTED;
    /** @type {string} */
    this.errors = '';

    for (const algorithm of ALGORITHMS) {
      this[algorithm] = Status.NOT_STARTED;
    }
  }

  setError(step, msg) {
    this[step] = Status.FAILURE;
    this.errors += msg + '\n';
  }
}

export class DatasetSummaries {
  constructor() {
    /** @type {DatasetSummary[]} */
    this.summaries = [];
    /** @type {Map<string, DatasetSummary>} */
    this.byId = new Map();
  }

  /**
   * Adds items at the end of the summaries
   *
   * @param  {...DatasetSummary} items
   */
  push(...items) {
    this.summaries.push(...items);
    items.forEach((item) => this.byId.set(item.id, item));
  }

  /**
   * Gets a dataset by id
   *
   * @param {string} id Dataset id
   * @returns The dataset
   */
  get(id) {
    return this.byId.get(id);
  }

  *items() {
    yield* this.summaries;
  }

  static async load(path) {
    const content = await readFile(path, { encoding: 'utf8' });
    const parsed = Papa.parse(content, {
      header: true,
      skipEmptyLines: 'greedy',
    });
    const items = parsed.data.map(DatasetSummary.fromRaw);
    const summaries = new DatasetSummaries();
    summaries.push(...items);
    return summaries;
  }

  async save(path) {
    const content = Papa.unparse({
      fields: ['id', 'downloaded', ...ALGORITHMS, 'errors'],
      data: this.summaries,
    });
    writeFile(path, content);
  }
}
