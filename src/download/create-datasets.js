import { Dataset } from '../dataset/dataset.js';
import { DatasetSummaries } from '../dataset/summary.js';
import { loadListing } from '../util/common.js';
import { Config } from '../util/config.js';
import { loadDatasetHandlers } from '../util/handler.js';
import { mapEntries } from '../util/iter.js';
import { setDownloaderRef, setSummaryRef } from './utils.js';

/**
 * Find a handler that supports the dataset
 *
 * @param {Dataset} dataset Dataset
 * @param {Iterable<[string, import('../util/handler.js').DatasetHandler]>} handlers Handler candidates
 * @returns Name of supporting handler
 */
function findHandler(dataset, handlers) {
  for (const [name, handler] of handlers) {
    if (handler.supports(dataset)) {
      return name;
    }
  }

  throw new Error('Unreachable - Default handler should have matched!');
}

/**
 * Initializes datasets
 *
 * @param {DatasetSummaries} summaries
 * @param {Config} config
 */
export async function createDatasets(summaries, config) {
  const listing = await loadListing(config);
  const handlers = await loadDatasetHandlers(config);
  const downloaders = new Map(mapEntries(handlers, (handler) => new handler.Downloader(config)));

  const datasets = listing.map((item) => new Dataset(item.id, config, item));
  for (const dataset of datasets) {
    const handler = findHandler(dataset, handlers);
    const summary = summaries.get(dataset.id);
    const downloader = downloaders.get(handler);

    dataset.handler = handler;
    setSummaryRef(dataset, summary);
    setDownloaderRef(dataset, downloader);
  }

  return datasets;
}
