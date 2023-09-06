export * from './downloader.js';

export function supports(dataset) {
  return /^hbm/i.test(dataset.id);
}
