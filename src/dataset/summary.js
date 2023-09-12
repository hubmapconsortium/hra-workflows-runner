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
  NOT_SUPPORTED: 'not supported',
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

  /**
   * Sets success for step and optionally clear errors
   *
   * @param {string} step
   * @param {boolean} [clearErrors=true]
   */
  setSuccess(step, clearErrors = true) {
    this[step] = Status.SUCCESS;
    if (clearErrors) {
      this.errors = '';
    }
  }

  /**
   * Sets failure for step and add an error message
   *
   * @param {string} step
   * @param {any} msg
   * @param {boolean} [mergeErrors=false]
   */
  setFailure(step, msg, mergeErrors = false) {
    const formattedMessage = JSON.stringify(msg)
      .slice(1, -1)
      .replace(/\\"/g, '"');
    const errors = mergeErrors ? this.errors + '\\n' : '';

    this[step] = Status.FAILURE;
    this.errors = errors + formattedMessage;
  }

  /**
   * Sets not supported for step
   *
   * @param {string} step
   */
  setNotSupported(step) {
    this[step] = Status.NOT_SUPPORTED;
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

  /**
   * Iterator of each summary item
   */
  *values() {
    yield* this.summaries;
  }

  /**
   * Load summaries from file
   *
   * @param {import('node:fs').PathLike} path File path
   */
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

  /**
   * Save summaries to file
   *
   * @param {import('node:fs').PathLike} path File path
   */
  async save(path) {
    const content = Papa.unparse({
      fields: ['id', 'downloaded', ...ALGORITHMS, 'errors'],
      data: this.summaries,
    });
    await writeFile(path, content);
  }
}
