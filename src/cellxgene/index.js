export * from './listing.js';
export * from './downloader.js';
export * from './job-generator.js';

export function supports(dataset) {
  return /cellxgene/i.test(dataset.id) && URL.canParse(dataset.id);
}
