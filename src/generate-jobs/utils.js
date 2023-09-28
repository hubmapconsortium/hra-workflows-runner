import { Dataset, createScratchGetSet } from '../dataset/dataset.js';
import { DatasetSummary } from '../dataset/summary.js';
import { getSummaryRef } from '../util/common.js';
import { IJobGenerator } from '../util/handler.js';

export const { get: getJobGeneratorRef, set: setJobGeneratorRef } =
  /** @type {import('../dataset/dataset.js').ScratchGetSetPair<IJobGenerator>} */
  (createScratchGetSet('jobGeneratorRef'));

/**
 * Get the summary instances for each datasets
 *
 * @param {Dataset[]} datasets
 * @returns {DatasetSummary[]}
 */
export function getSummaryRefMany(datasets) {
  return datasets.map(getSummaryRef).filter((ref) => ref !== undefined);
}
