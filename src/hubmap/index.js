export * from './listing.js';
export * from './downloader.js';
export * from '../xconsortia/job-generator.js';

export function supports(dataset) {
  return /^hbm/i.test(dataset.id);
}
