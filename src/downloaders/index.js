import CellXGeneDownloader from './cellxgene.js';
import Downloader from './downloader.js';
import GtexDownloader from './gtex.js';
import HubmapDownloader from './hubmap.js';

/**
 * Available downloaders
 *
 * @type {(typeof Downloader)[]}
 */
const DOWNLOADER_CLASSES = [
  CellXGeneDownloader,
  GtexDownloader,
  HubmapDownloader,
];

/**
 * Creates a downloader for a dataset
 * @param {object} dataset The dataset
 * @param {string} outDir Output directory
 * @param {string} cacheDir Temporary cache directory
 * @returns {Downloader | undefined} The downloader or undefined if the dataset is not supported
 */
export function createDownloader(dataset, outDir, cacheDir) {
  const cls = DOWNLOADER_CLASSES.find((cls) => cls.supports(dataset));
  return cls && new cls(dataset, outDir, cacheDir);
}
