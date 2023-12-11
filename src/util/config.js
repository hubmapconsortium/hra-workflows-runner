import { env } from 'node:process';

/** Indicator that a configuration value is required */
export const REQUIRED = Symbol('Required config value');

/** Manage global/environment configuration */
export class Config {
  constructor() {
    /** @type {Map<string, any>} */
    this.config = new Map();
  }

  /**
   * Fetch a configuration value
   *
   * @template T
   * @param {string} key Configuration key
   * @param {T} [defaultValue] Default value
   * @returns {T | undefined}
   * @throws If the key doesn't exist and no default value was provided
   */
  get(key, defaultValue) {
    // Checking arguments.length allows undefined to be used as a default value
    if (arguments.length < 2) {
      defaultValue = REQUIRED;
    }

    if (this.config.has(key)) {
      return this.config.get(key);
    } else if (defaultValue !== REQUIRED) {
      return defaultValue;
    } else {
      throw new Error(`Missing required config value for '${key}'`);
    }
  }

  /**
   * Associate a new value with a configuration key
   *
   * @template T
   * @param {string} key Configuration key
   * @param {T} value New value
   */
  set(key, value) {
    this.config.set(key, value);
  }

  /**
   * Gets all keys that starts with a prefix
   *
   * @param {string} prefix Key prefix
   */
  getPrefixedKeys(prefix) {
    const keys = Array.from(this.config.keys());
    return keys.filter((key) => key.startsWith(prefix));
  }

  /**
   * Checks whether all required keys exist in the config
   *
   * @param {string[]} required Required keys
   * @throws If not all keys are present
   */
  validate(required) {
    for (const key of required) {
      if (!this.config.has(key)) {
        throw new Error(`Missing configuration key '${key}'`);
      }
    }

    return this;
  }

  /**
   * Loads the environment into this configuration object
   *
   * @param {function(string, string): any} [reviver] Reviver function
   */
  loadEnv(reviver = (_, val) => val) {
    for (const [key, value] of Object.entries(env)) {
      this.set(key, reviver(key, value));
    }

    return this;
  }
}
