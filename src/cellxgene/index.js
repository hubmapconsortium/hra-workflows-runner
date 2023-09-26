export * from './downloader.js';

export function supports(dataset) {
  return /cellxgene/i.test(dataset.id) && URL.canParse(dataset.id);
}
