import { readFile, writeFile } from 'node:fs/promises';
import { Config } from '../util/config.js';
import {
  getDataDir,
  getDataFilePath,
  getDatasetFilePath,
  getDirForId,
  getJobFilePath,
  getJobFilePathWithSuffix,
} from '../util/paths.js';

/** Non-serialized config property */
const CONFIG = Symbol('Configuration reference');
/** Non-serialized scratch data property */
const SCRATCH = Symbol('Non-serialized data');

/**
 * @template T
 * @typedef ScratchGetSetPair
 * @property {function(Dataset): T} get Getter function
 * @property {function(Dataset, T): T} set Setter function
 */

/**
 * Creates a get/set pair for a scratch key
 *
 * @template T
 * @param {string|number|symbol} key Scratch key
 * @returns {ScratchGetSetPair<T>}
 */
export function createScratchGetSet(key) {
  return {
    get: (dataset) => dataset.scratch[key],
    set: (dataset, value) => (dataset.scratch[key] = value),
  };
}

/** Dataset metadata class */
export class Dataset {
  constructor(id, config, extra = {}) {
    /** @type {string} */
    this.id = id;
    /** @type {string} */
    this.handler = '';
    /** @type {string} */
    this.organ = '';
    /** @type {Config} */
    this[CONFIG] = config;
    /** @type {Record<string, any>} */
    this[SCRATCH] = {};

    Object.assign(this, extra);
  }

  /** Reference to the global configuration */
  get config() {
    return this[CONFIG];
  }

  /** Reference to non-serialized scratch object */
  get scratch() {
    return this[SCRATCH];
  }

  /** Directory name for this dataset */
  get dir() {
    return getDirForId(this.id);
  }

  /** Directory path for this dataset */
  get dirPath() {
    return getDataDir(this.config, this.dir);
  }

  /** File path to dataset file */
  get filePath() {
    return getDatasetFilePath(this.config, this.dir);
  }

  /** File path to data file */
  get dataFilePath() {
    return getDataFilePath(this.config, this.dir);
  }

  /** File path to job file */
  get jobFilePath() {
    return getJobFilePath(this.config, this.dir);
  }

  /**
   * Get the job file path with an added suffix
   *
   * @param {string} suffix Suffix
   */
  jobFilePathWithSuffix(suffix) {
    return getJobFilePathWithSuffix(this.config, this.dir, suffix);
  }

  /**
   * Load dataset from file
   *
   * @param {import('node:fs').PathLike} path File path
   * @param {Config} config Configuration
   */
  static async load(path, config) {
    const content = await readFile(path, { encoding: 'utf8' });
    const data = JSON.parse(content);
    return new Dataset(data.id, config, data);
  }

  /**
   * Saves dataset to file
   *
   * @param {import('node:fs').PathLike} [path] File path
   */
  async save(path = this.filePath) {
    const content = JSON.stringify(this, undefined, 2);
    await writeFile(path, content, { encoding: 'utf8' });
  }
}
