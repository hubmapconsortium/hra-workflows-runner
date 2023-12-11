/** Used to give each timer a unique id */
let timerCounter = 0;

/**
 * Logs the start, end (or error), and execution time of an usually async callback.
 * 
 * @template R Result type
 * @param {string} event Event name
 * @param  {[...any, () => R | Promise<R>]} args Event arguments and callback
 * @returns {Promise<R>} The value returned by the callback
 */
export async function logEvent(event, ...args) {
  const [exec] = args.splice(-1, 1);
  const timerId = `${event}:Timer(${timerCounter++}) -- Args: ${args.join(' ')}`;
  try {
    console.debug(`${event}:Start -- Args:`, ...args);
    console.time(timerId);
    const result = await exec();
    console.debug(`${event}:End -- Args:`, ...args);
    return result;
  } catch (error) {
    console.error(`${event}:Failure -- Args:`, ...args, `-- Cause: ${error.message}`);
    throw error;
  } finally {
    console.timeEnd(timerId);
  }
}
