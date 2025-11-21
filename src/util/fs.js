import { access, constants, mkdir, open, writeFile } from 'node:fs/promises';
import { install, Agent } from 'undici';
import { concurrentMap } from './concurrent-map.js';

install();

/**
 * Ensure that all specified directories are created
 *
 * @param  {...import('node:fs').PathLike} dirs Directory paths
 */
export async function ensureDirsExist(...dirs) {
  const createDir = async (dir) => await mkdir(dir, { recursive: true });
  await concurrentMap(dirs, createDir);
}

/**
 * Tests whether a file exists
 *
 * @param {import('node:fs').PathLike} path File path
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
 * Check whether a fetch response was successful otherwise throw an error
 *
 * @param {Response} response The response object
 * @param {string} [msg] Message added before the status and text
 */
export function checkFetchResponse(response, msg = 'Download failed') {
  if (!response.ok) {
    const { status, statusText } = response;
    throw new Error(`${msg}: ${status}:${statusText}`);
  }
}

/**
 * Download a remote file
 *
 * @param {import('node:fs').PathLike} dest Destination file
 * @param {string | URL} src Source url
 * @param {RequestInfo & { overwrite?: boolean }} [options] Additional options
 */
export async function downloadFile(dest, src, options = {}) {
  const fileHandle = await openWriteFile(dest, options);
  if (fileHandle === undefined) {
    return;
  }

  try {
    const resp = await fetch(src, {
      dispatcher: new Agent({ connectTimeout: 300e3 }),
      ...options,
    });
    checkFetchResponse(resp);

    await writeFile(fileHandle, resp.body, { encoding: 'utf8' });
  } finally {
    await fileHandle.close();
  }
}

/**
 * Opens a file for writing
 *
 * @param {string} path Path to file
 * @param {{ overwrite?: boolean }} options Options controlling how to open the file
 * @returns A file handle or undefined if overwrite is false and the file already exists
 */
async function openWriteFile(path, options) {
  const { overwrite = false } = options;
  const flags = overwrite ? 'w' : 'wx';

  try {
    return await open(path, flags);
  } catch (error) {
    if (!overwrite && error.code === 'EEXIST') {
      return undefined;
    }

    throw error;
  }
}
