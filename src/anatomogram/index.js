export * from './listing.js';
export * from './downloader.js';
export * from './job-generator.js';

export function supports(dataset) {
  return /^anatomogram/i.test(dataset.id);
}
