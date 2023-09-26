/**
 * @template T, V
 * @extends Map<T, V>
 */
export class Cache extends Map {
  /**
   * Gets a value or set a default
   *
   * @param {K} key Key
   * @param {V} defaultValue Default value if key is not set
   * @returns {V}
   */
  setDefault(key, defaultValue) {
    if (!this.has(key)) {
      this.set(key, defaultValue);
    }

    return this.get(key);
  }

  /**
   * Gets a value or sets a default produced by the provided function
   *
   * @param {K} key Key
   * @param {function(): V} fn Function to call if key is not set
   * @returns {V}
   */
  setDefaultFn(key, fn) {
    if (!this.has(key)) {
      this.set(key, fn());
    }

    return this.get(key);
  }
}
