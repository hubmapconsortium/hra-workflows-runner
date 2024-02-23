import { join } from 'node:path';
import { loadJson } from '../util/common.js';
import { concurrentMap } from '../util/concurrent-map.js';
import { Config } from '../util/config.js';
import { ALGORITHMS, DEFAULT_MAX_CONCURRENCY, FORCE, MAX_CONCURRENCY } from '../util/constants.js';
import { downloadFile, ensureDirsExist } from '../util/fs.js';
import { getOutputDir } from '../util/paths.js';

/**
 * @typedef {Record<string, string | Record<string, any>>} RawOrganMetadata
 */

/** Template for organ metadata file urls */
const ORGAN_METADATA_URL_TEMPLATE =
  'https://raw.githubusercontent.com/hubmapconsortium/hra-workflows/main/containers/{{algorithm}}/context/organ-metadata.json';

/**
 * Raw metadata loaded from file cached by algorithm
 * @type {Map<string, Promise<RawOrganMetadata>>}
 */
const cachedMetadataFileDownload = new Map();

/**
 * Tests whether a value is a string
 *
 * @param {unknown} value Value to test
 * @returns {value is string}
 */
function isString(value) {
  return typeof value === 'string';
}

/**
 * Tests whether a value is an object
 *
 * @param {unknown} value Value to test
 * @returns {value is Record<string, any>}
 */
function isObject(value) {
  return typeof value === 'object';
}

export class OrganMetadata {
  /**
   * Loads organ metadata for a specific algorithm
   *
   * @param {string} algorithm Algorithm for which to load metadata
   * @param {Config} config Configuration
   */
  static async load(algorithm, config) {
    if (!cachedMetadataFileDownload.has(algorithm)) {
      const url = ORGAN_METADATA_URL_TEMPLATE.replace('{{algorithm}}', algorithm);
      const dir = join(getOutputDir(config), 'organ-metadata');
      const file = join(dir, `${algorithm}.json`);
      const download = async () => {
        await ensureDirsExist(dir);
        await downloadFile(file, url, { overwrite: config.get(FORCE, false) });
        return await loadJson(file);
      };

      cachedMetadataFileDownload.set(algorithm, download());
    }

    const result = await cachedMetadataFileDownload.get(algorithm);
    return new OrganMetadata({ ...result });
  }

  constructor(metadata) {
    /** @type {RawOrganMetadata} */
    this.metadata = metadata;
  }

  /**
   * All organ codes in the metadata
   */
  get organs() {
    return Object.keys(this.metadata);
  }

  /**
   * Gets the metadata for an organ.
   * May return a redirect organ code.
   *
   * @param {string} organ Organ code
   * @returns {RawOrganMetadata[string] | undefined}
   */
  get(organ) {
    return this.metadata[organ];
  }

  /**
   * Test whether an organ is in the metadata
   *
   * @param {string} organ Organ code
   */
  has(organ) {
    return this.get(organ) !== undefined;
  }

  /**
   * Resolves an organ following all redirect organ codes
   *
   * @param {string} organ Organ code
   * @returns {string} The resolved organ code
   */
  resolve(organ) {
    const value = this.get(organ);
    return isString(value) ? this.resolve(value) : organ;
  }

  /**
   * Computes the specificity of an organ.
   * The specificity of an organ is the number of redirect organ codes
   * that has to be followed to resolve the organ.
   *
   * @param {string} organ Organ code
   * @param {number} initial Initial specificity value
   * @returns {number}
   */
  getSpecificity(organ, initial = 0) {
    const value = this.get(organ);
    return isString(value) ? this.getSpecificity(value, initial + 1) : initial;
  }
}

export class OrganMetadataCollection {
  /**
   * Loads all organ metadata for every algorithm
   *
   * @param {Config} config Configuration
   */
  static async load(config) {
    const collection = await concurrentMap(ALGORITHMS, (algorithm) => OrganMetadata.load(algorithm, config), {
      maxConcurrency: config.get(MAX_CONCURRENCY, DEFAULT_MAX_CONCURRENCY),
    });

    return new OrganMetadataCollection(collection);
  }

  constructor(collection) {
    /** @type {OrganMetadata[]} */
    this.collection = collection;
  }

  /**
   * All organ codes in the metadata
   */
  get organs() {
    return Array.from(new Set(this.collection.flatMap((meta) => meta.organs)));
  }

  /**
   * Gets the metadata for an organ.
   * May return a redirect organ code.
   *
   * @param {string} organ Organ code
   */
  get(organ) {
    return this.collection.map((meta) => meta.get(organ));
  }

  /**
   * Test whether an organ is in the metadata
   *
   * @param {string} organ Organ code
   */
  has(organ) {
    return this.collection.some((meta) => meta.has(organ));
  }

  /**
   * Resolves an organ following redirect organ codes.
   * Resolving stops when at least one of the organ metadata
   * in the collection can not resolve further.
   *
   * @param {string} organ Organ code
   * @returns {string} The resolved organ code
   */
  resolve(organ) {
    const values = this.get(organ);
    if (values.some(isObject)) {
      return organ;
    }

    const next = values.find(isString);
    return next !== undefined ? this.resolve(next) : organ;
  }

  /**
   * Computes the specificity of an organ.
   * The specificity of an organ is the number of redirect organ codes
   * that has to be followed to resolve the organ.
   *
   * @param {string} organ Organ code
   * @param {number} initial Initial specificity value
   * @returns The highest specificity reported by the metadata in the collection
   */
  getSpecificity(organ, initial = 0) {
    const specificities = this.collection.map((meta) => meta.getSpecificity(organ, initial));
    return Math.max(...specificities);
  }

  /**
   * Selects the organ with the highest specificity
   *
   * @param  {...string} organs Organ codes to select among
   */
  selectBySpecificity(...organs) {
    organs.sort((a, b) => this.getSpecificity(b) - this.getSpecificity(a));
    return organs[0];
  }
}
