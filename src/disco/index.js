export * from './listing.js';
export * from './downloader.js';
export * from './job-generator.js';

export function supports(dataset) {
  return /^disco/i.test(dataset.id);
}
