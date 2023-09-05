import { env } from 'node:process';

export const REQUIRED = Symbol('Required config value');

export class Config {
  constructor() {
    /** @type {Map<string, any>} */
    this.config = new Map();
  }

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

  set(key, value) {
    this.config.set(key, value);
  }

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

  loadEnv(reviver = (_, val) => val) {
    for (const [key, value] of Object.entries(env)) {
      this.set(key, reviver(key, value));
    }

    return this;
  }
}
