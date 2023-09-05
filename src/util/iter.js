/**
 * Maps an iterable of key-value pairs
 *
 * @template K
 * @template V
 * @template V2
 * @param {Iterable<[K, V]>} iterable Iterable object
 * @param {function(V, K): V2} callbackFn Mapping function
 */
export function* mapEntries(iterable, callbackFn) {
  for (const [key, value] of iterable) {
    yield /** @type {[K, V2]} */ ([key, callbackFn(value, key)]);
  }
}

/**
 * Splits an iterable into groups
 *
 * @template K
 * @template V
 * @param {Iterable<V>} iterable Iterable object
 * @param {function(V): K} keyFn Key function
 */
export function groupBy(iterable, keyFn) {
  /** @type {Map<K, V[]>} */
  const groups = new Map();
  for (const value of iterable) {
    const key = keyFn(value);
    const group = setDefault(groups, key, []);
    group.push(value);
  }

  return groups;
}

/**
 * Gets the current value for a key.
 * If it doesn't exist set and returns a default value.
 * 
 * @template K
 * @template V
 * @param {Map<K, V>} map 
 * @param {K} key 
 * @param {V} defaultValue 
 * @returns {V}
 */
function setDefault(map, key, defaultValue) {
  if (map.has(key)) {
    return map.get(key);
  }

  map.set(key, defaultValue);
  return defaultValue;
}

/**
 * Simple diff. Removes all values in an iterable that are also in the exclude iterable
 * 
 * @template T
 * @param {Iterable<T>} iterable Iterable object
 * @param {Iterable<T>} exclude Excluded values
 */
export function* diff(iterable, exclude) {
  const excludeSet = new Set(exclude);
  for (const value of iterable) {
    if (!excludeSet.has(value)) {
      yield value;
    }
  }
}
