import { access, constants, mkdir, open, writeFile } from 'node:fs/promises';
import { concurrentMap } from './concurrent-map.js';

export async function ensureDirsExist(...dirs) {
  const createDir = async (dir) => await mkdir(dir, { recursive: true });
  await concurrentMap(dirs, createDir);
}

export async function fileExists(path) {
  try {
    await access(path, constants.R_OK | constants.W_OK);
    return true;
  } catch {
    return false;
  }
}

export async function downloadFile(dest, src, options) {
  const fileHandle = await openWriteFile(dest, options);
  if (fileHandle === undefined) {
    return;
  }

  try {
    const resp = await fetch(src);
    if (!resp.ok) {
      const { status, statusText } = resp;
      const message = `Download failed: ${status}:${statusText}`;
      throw new Error(message);
    }

    await writeFile(fileHandle, resp.body, { encoding: 'utf8' });
  } finally {
    await fileHandle.close();
  }
}

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
