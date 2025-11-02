import { createReadStream } from 'fs';
import { createGunzip } from 'zlib';
import Papa from 'papaparse';

export function normalizeCsvUrl(url) {
  if (url.startsWith('https://docs.google.com/spreadsheets/d/') && url.indexOf('export?format=csv') === -1) {
    const splitUrl = url.split('/');
    if (splitUrl.length === 7) {
      const sheetId = splitUrl[5];
      const gid = splitUrl[6].split('=')[1];
      return `https://docs.google.com/spreadsheets/d/${sheetId}/export?format=csv&gid=${gid}`;
    }
  }
  return url;
}

export async function* readLines(inputFile) {
  let inputStream = !inputFile || inputFile === '-' ? process.stdin : createReadStream(inputFile, { autoClose: true });
  if (inputFile?.endsWith('.gz')) {
    inputStream = inputStream.pipe(createGunzip());
  }
  const decoder = new TextDecoder('utf-8');
  let buffer = '';
  for await (const chunk of inputStream) {
    buffer += decoder.decode(chunk, { stream: true });
    let lines = buffer.split('\n');
    buffer = lines.pop();
    for (const line of lines) {
      yield line;
    }
  }
  if (buffer.length > 0) {
    yield buffer;
  }
}

export async function* readCsv(input, options = { skipEmptyLines: true, header: true }) {
  if (options.header) {
    let header;
    const newOpts = { ...options, header: false };
    for await (const line of readLines(input)) {
      const row = Papa.parse(line, newOpts)?.data ?? [];
      if (row) {
        if (!header) {
          header = row[0];
        } else {
          const result = {};
          for (let i = 0; i < header.length; i++) {
            result[header[i]] = row[0][i];
          }
          yield result;
        }
      }
    }
  } else {
    return createReadStream(input).pipe(Papa.parse(Papa.NODE_STREAM_INPUT, options));
  }
}
