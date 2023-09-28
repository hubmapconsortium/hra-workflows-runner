import { readFile, writeFile } from 'node:fs/promises';
import { Config } from '../util/config.js';
import {
  getDataDir,
  getDataFilePath,
  getDatasetFilePath,
  getDirForId,
  getJobFilePath,
} from '../util/paths.js';

const CONFIG = Symbol('Configuration reference');
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

  get config() {
    return this[CONFIG];
  }

  get scratch() {
    return this[SCRATCH];
  }

  get dir() {
    return getDirForId(this.id);
  }

  get dirPath() {
    return getDataDir(this.config, this.dir);
  }

  get filePath() {
    return getDatasetFilePath(this.config, this.dir);
  }

  get dataFilePath() {
    return getDataFilePath(this.config, this.dir);
  }

  get jobFilePath() {
    return getJobFilePath(this.config, this.dir);
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
