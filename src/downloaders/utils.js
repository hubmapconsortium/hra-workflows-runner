import { createWriteStream } from 'node:fs';
import { access, constants } from 'node:fs/promises';
import { IncomingMessage } from 'node:http';
import { get } from 'node:https';
import { Writable } from 'node:stream';
import { pipeline } from 'node:stream/promises';

// Unique object
const NULL = {};

/**
 * Wraps a function to only run once.
 * All subsequent calls return the same result as the first call.
 * @template {function} T
 * @param {T} fun Function to wrap
 * @returns {T} The wrapped function
 */
export function once(fun) {
  let result = NULL;
  return function (...args) {
    return result !== NULL ? result : (result = fun.apply(this, args));
  };
}

/**
 * Test whether a file exists
 * @param {string} path File path to test
 * @returns {Promise<boolean>} A promise that resolves to true if the file exists
 */
export async function fileExists(path) {
  try {
    await access(path, constants.R_OK | constants.W_OK);
    return true;
  } catch {
    return false;
  }
}

/**
 * @typedef FileWriteStreamOptions
 * @type {object}
 * @property {boolean} [overwrite=false] Whether to overwrite an existing file
 */

/**
 * Open a stream for writing to a file
 * @param {string} path File path
 * @param {FileWriteStreamOptions} [options] Additional options
 * @returns {Promise<Writable | undefined>} A promise that resolves to the stream once it is ready for writes
 */
export async function openFileWriteStream(path, options = {}) {
  const { overwrite = false } = options;
  const flags = overwrite ? 'w' : 'wx';
  return new Promise((resolve, reject) => {
    const stream = createWriteStream(path, { flags })
      .once('ready', () => resolve(stream))
      .once('error', (error) =>
        !overwrite && error.code === 'EEXIST'
          ? resolve(undefined)
          : reject(error)
      );
  });
}

/**
 * Open a read stream to a remote resource
 * @param {string} url Url of resource to fetch
 * @returns {Promise<IncomingMessage>} A promise that resolves to the stream
 */
export async function openUrlReadStream(url) {
  return new Promise((resolve, reject) => {
    get(url, (incoming) => {
      if (incoming.statusCode !== 200) {
        const { statusCode, statusMessage } = incoming;
        const message = `Download failed: ${statusCode}:${statusMessage} from ${url}`;
        reject(new Error(message));
      } else {
        resolve(incoming);
      }
    }).once('error', reject);
  });
}

/**
 * Download a remote file
 * @param {string} url Url of remote file
 * @param {string} localPath Local file path
 * @param {FileWriteStreamOptions} [options] Additional options
 */
export async function downloadRemoteFile(url, localPath, options = {}) {
  const fileStream = await openFileWriteStream(localPath, options);
  if (fileStream) {
    const urlStream = await openUrlReadStream(url);
    await pipeline(urlStream, fileStream);
  }
}
