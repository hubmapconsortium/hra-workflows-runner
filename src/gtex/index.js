export * from './downloader.js';

export function supports(dataset) {
  return /^gtex/i.test(dataset.id);
}
