import { writeFile } from 'node:fs/promises';
import Papa from 'papaparse';
import { loadCsv } from '../util/common.js';
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

/**
 * Builtin steps
 * @readonly
 * @enum {string}
 */
export const Step = {
  DOWNLOADED: 'downloaded',
};

/**
 * Builtin algorithm steps
 * @readonly
 * @enum {string}
 */
export const AlgorithmStep = ALGORITHMS.reduce(
  (result, algorithm) => ({
    ...result,
    [algorithm.toUpperCase()]: algorithm,
  }),
  /** @type {Object<string, string>} */ ({})
);

const STEPS = Object.values({
  ...Step,
  ...AlgorithmStep,
});

export class DatasetSummary {
  /** @private */
  static fromRaw(data) {
    return Object.assign(new DatasetSummary(data.id), data);
  }

  /**
   * Get the property storing errors for a step
   *
   * @param {string} step
   */
  static getErrorPropertyForStep(step) {
    return `errors:${step}`;
  }

  /**
   * Creates a new dataset summary item
   *
   * @param {string} id Dataset id
   */
  constructor(id) {
    /** @type {string} */
    this.id = id;
    /** @type {string} */
    this.info = ''

    for (const step of STEPS) {
      this.setStatus(step, Status.NOT_STARTED);
      this.setError(step, '');
    }
  }

  /**
   * Get the status of a step
   *
   * @param {string} step
   * @returns {Status}
   */
  getStatus(step) {
    return this[step] ?? Status.NOT_STARTED;
  }

  /**
   * Updates the status of a step
   *
   * @param {string} step
   * @param {Status} status
   */
  setStatus(step, status) {
    this[step] = status;
  }

  /**
   * Get any errors for a step
   *
   * @param {string} step
   * @returns {string}
   */
  getError(step) {
    const prop = DatasetSummary.getErrorPropertyForStep(step);
    return this[prop];
  }

  /**
   * Sets the error message for a step
   *
   * @param {string} step
   * @param {any} message Error message
   */
  setError(step, message) {
    const prop = DatasetSummary.getErrorPropertyForStep(step);
    const formattedMessage = JSON.stringify(message)
      .slice(1, -1)
      .replace(/\\"/g, '"');

    this[prop] = formattedMessage;
  }

  /**
   * Sets success for step and optionally clear errors
   *
   * @param {string} step
   */
  setSuccess(step) {
    this.setStatus(step, Status.SUCCESS);
    this.setError(step, '');
  }

  /**
   * Sets success for step on multiple summaries
   *
   * @param {DatasetSummary[]} summaries
   * @param {string} step
   */
  static setSuccessMany(summaries, step) {
    summaries.forEach((summary) => summary.setSuccess(step));
  }

  /**
   * Sets failure for step and add an error message
   *
   * @param {string} step
   * @param {any} message
   */
  setFailure(step, message) {
    this.setStatus(step, Status.FAILURE);
    this.setError(step, message);
  }

  /**
   * Sets failure for step on multiple summaries
   *
   * @param {DatasetSummary[]} summaries
   * @param {string} step
   * @param {any} msg
   * @param {boolean} [mergeErrors=false]
   */
  static setFailureMany(summaries, step, msg, mergeErrors) {
    summaries.forEach((summary) => summary.setFailure(step, msg, mergeErrors));
  }

  /**
   * Sets not supported for step
   *
   * @param {string} step
   */
  setNotSupported(step) {
    this.setStatus(step, Status.NOT_SUPPORTED);
    this.setError(step, '');
  }

  /**
   * Sets not supported for step on multiple summaries
   *
   * @param {DatasetSummary[]} summaries
   * @param {string} step
   */
  static setNotSupportedMany(summaries, step) {
    summaries.forEach((summary) => summary.setNotSupported(step));
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
   * Filters summary items whose step has the specified status
   *
   * @param {string} step
   * @param {Status} status
   */
  filterByStatus(step, status) {
    return this.summaries.filter(
      (summary) => summary.getStatus(step) === status
    );
  }

  /**
   * Load summaries from file
   *
   * @param {import('node:fs').PathLike} path File path
   */
  static async load(path) {
    const rows = await loadCsv(path);
    const items = rows.map(DatasetSummary.fromRaw);
    const summaries = new DatasetSummaries();
    summaries.push(...items);
    return summaries;
  }

  /**
   * Save summaries to file
   *
   * @param {import('node:fs').PathLike} path File path
   * @param {string[]} [steps] Steps to serialize
   */
  async save(path, steps = STEPS) {
    const errorProps = steps.map(DatasetSummary.getErrorPropertyForStep);
    const content = Papa.unparse({
      fields: ['id', ...steps, ...errorProps, 'info'],
      data: this.summaries,
    });
    await writeFile(path, content);
  }
}
