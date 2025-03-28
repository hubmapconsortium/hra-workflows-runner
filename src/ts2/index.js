export * from './listing.js';
export * from './downloader.js';
export * from './job-generator.js';

export function supports(dataset) {
  return /^ts2/i.test(dataset.id);
}
