import { join } from 'node:path';
import { concurrentMap } from '../util/concurrent-map.js';
import { FORCE } from '../util/constants.js';
import { checkFetchResponse, downloadFile, ensureDirsExist } from '../util/fs.js';
import { getCacheDir } from '../util/paths.js';

export async function getCollections(figshareId) {
  const url = new URL(`https://api.figshare.com/v2/articles/${figshareId}/files?limit=100`);
  const resp = await fetch(url, { method: 'GET' });
  checkFetchResponse(resp, 'TS2: Failed to fetch list of collections');

  const collections = await resp.json();
  return collections;
}

export async function cacheCollections(figshareId, config) {
  const collections = await getCollections(figshareId);
  const dataDir = join(getCacheDir(config), 'ts2');
  await ensureDirsExist(dataDir);

  await concurrentMap(collections, (collection) => {
    const dataFilePath = join(dataDir, collection.name);
    const downloadUrl = collection.download_url;
    return downloadFile(dataFilePath, downloadUrl, {
      overwrite: config.get(FORCE, false),
    });
  });

  return collections;
}
