/**
 * Maps each value in an array concurrently
 *
 * @template T Value type
 * @template R Mapped value type
 * @param {T[]} array Original array
 * @param {function(T, number, T[]): R | PromiseLike<R>} mapper Potentially asynchronous mapping function
 * @param {{maxConcurrency?: number}} [options] Additional options
 * @returns {Promise<R[]>} The mapped array
 */
export async function concurrentMap(array, mapper, options = {}) {
  const executor = new ConcurrentMapExecutor(array, mapper, options);
  await executor.run();
  return executor.result;
}

/**
 * Class implementing the concurrentMap logic
 *
 * @template T Input value type
 * @template R Result value type
 */
class ConcurrentMapExecutor {
  constructor(array, mapper, options) {
    /** @type {T[]} */
    this.array = array;
    /** @type {function(T, number, T[]): R | PromiseLike<R>} */
    this.mapper = mapper;
    /** @type {{maxConcurrency?: number}} */
    this.options = options;
    /** @type {number} */
    this.index = 0;
    /** @type {boolean} */
    this.errored = false;
    /** @type {R[]} */
    this.result = Array(array.length).fill(undefined);
  }

  /**
   * Asynchrously execute the mapper function over all the values
   */
  async run() {
    const {
      array: { length },
      options: { maxConcurrency = 10 },
    } = this;
    const numExecutors = Math.min(length, maxConcurrency);
    const executors = Array(numExecutors)
      .fill(undefined)
      .map(() => this.execute());

    await Promise.all(executors);
  }

  /**
   * Maps a single value at a time until all values have been mapped or an error is thrown
   */
  async execute() {
    const { array, mapper, result } = this;
    while (this.index < array.length && !this.errored) {
      const index = this.index++;
      const item = array[index];
      try {
        result[index] = await mapper(item, index, array);
      } catch (error) {
        this.errored = true;
        throw error;
      }
    }
  }
}
