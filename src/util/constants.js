// Environment variables
export const FORCE = 'FORCE';
export const MAX_CONCURRENCY = 'MAX_CONCURRENCY';
export const DATASET = 'DATASET';
export const VERSION = 'VERSION';
export const OUTPUT_DIR = 'OUTPUT_DIR';
export const DATA_REPO_DIR = 'DATA_REPO_DIR';
export const CACHE_DIR = 'CACHE_DIR';
export const SRC_DIR = 'SRC_DIR';
export const DATASET_LIST_URL = 'DATASET_LIST_URL';
export const DATASET_ID_COLUMN = 'DATASET_ID_COLUMN';
export const DATASET_LINK_COLUMN = 'DATASET_LINK_COLUMN';

export const REQUIRED_ENV_VARIABLES = [
  OUTPUT_DIR,
  DATA_REPO_DIR,
  SRC_DIR,
  DATASET_LIST_URL,
  DATASET_ID_COLUMN,
  DATASET_LINK_COLUMN,
];

// File names
export const LISTING_FILE = 'listing.csv';
export const SUMMARIES_FILE = 'summaries.csv';
export const DATA_FILE = 'data.h5ad';
export const ALGORITHM_REPORT_FILE = 'report.json';

// Algorithms
export const ALGORITHMS = ['azimuth', 'celltypist', 'popv'];

// Default values
export const DEFAULT_MAX_CONCURRENCY = 10;
export const DEFAULT_CACHE_DIR = './tmp';
