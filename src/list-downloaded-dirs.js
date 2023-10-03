import { Status, Step } from './dataset/summary.js';
import { getConfig, loadSummaries } from './util/common.js';
import { getDataDir, getDirForId } from './util/paths.js';

async function main() {
  const config = getConfig();
  const summaries = await loadSummaries(config);
  const downloadedItems = summaries.filterByStatus(
    Step.DOWNLOADED,
    Status.SUCCESS
  );

  for (const item of downloadedItems) {
    const dir = getDataDir(config, getDirForId(item.id));
    console.log(dir);
  }
}

main();
