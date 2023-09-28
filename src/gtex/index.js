export * from './downloader.js';
export * from './job-generator.js';

export function supports(dataset) {
  return /^gtex/i.test(dataset.id);
}
