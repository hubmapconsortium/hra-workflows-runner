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

class ConcurrentMapExecutor {
  constructor(array, mapper, options) {
    this.array = array;
    this.mapper = mapper;
    this.options = options;
    this.index = 0;
    this.errored = false;
    this.result = Array(array.length).fill(undefined);
  }

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
