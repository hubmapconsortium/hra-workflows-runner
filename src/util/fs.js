import { open, writeFile } from 'node:fs/promises';

export async function downloadFile(dest, src, options) {
  const fileHandle = await openWriteFile(dest, options);
  if (fileHandle === undefined) {
    return;
  }

  const resp = await fetch(src);
  if (!resp.ok) {
    const { status, statusText } = resp;
    const message = `Download failed: ${status}:${statusText}`;
    throw new Error(message);
  }

  const data = await resp.text();
  await writeFile(fileHandle, data, { encoding: 'utf8' });
}

async function openWriteFile(path, options) {
  const { overwrite = false } = options;
  const flags = overwrite ? 'w' : 'wx';

  try {
    return open(path, flags);
  } catch (error) {
    if (!overwrite && error.code === 'EEXIST') {
      return undefined;
    }

    throw error;
  }
}
