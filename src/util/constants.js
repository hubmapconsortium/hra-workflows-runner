// Environment variables
export const FORCE = 'FORCE';
export const MAX_CONCURRENCY = 'MAX_CONCURRENCY';
export const PYTHON_LOG_LEVEL = 'PYTHON_LOG_LEVEL';
export const DATASET_HANDLERS = 'DATASET_HANDLERS';
export const DATASET = 'DATASET';
export const VERSION = 'VERSION';
export const DATASETS_DIR = 'DATASETS_DIR';
export const CROSSWALKING_TABLES_DIR = 'CROSSWALKING_TABLES_DIR';
export const OUTPUT_DIR = 'OUTPUT_DIR';
export const DATA_REPO_DIR = 'DATA_REPO_DIR';
export const CACHE_DIR = 'CACHE_DIR';
export const MODELS_DIR = 'MODELS_DIR';
export const SRC_DIR = 'SRC_DIR';
export const DATASET_LIST = 'DATASET_LIST';
export const DATASET_LIST_URL = 'DATASET_LIST_URL';
export const DATASET_COLUMN_ID = 'DATASET_COLUMN_ID';
export const DATASET_MIN_CELL_COUNT = 'DATASET_MIN_CELL_COUNT';
export const DATASET_DATA_MAX_SIZE = 'DATASET_DATA_MAX_SIZE';
export const POPV_METHOD_COUNT = 'POPV_METHOD_COUNT';

export const REQUIRED_ENV_VARIABLES = [DATASET, DATASETS_DIR, OUTPUT_DIR, DATA_REPO_DIR, MODELS_DIR, SRC_DIR];

// File names
export const LISTING_FILE = 'listing.csv';
export const SUMMARIES_FILE = 'summaries.csv';
export const DATASET_FILE = 'dataset.json';
export const DATA_FILE = 'data.h5ad';
export const JOB_FILE = 'job.json';
export const ALGORITHM_REPORT_FILE = 'report.json';
export const ALGORITHM_SUMMARY_JSON_LD_FILE = 'summary.jsonld';
export const ALGORITHM_ANNOTATIONS_FILE = 'annotations.csv';

// Algorithms
export const ALGORITHMS = ['azimuth', 'celltypist', 'popv', 'frmatch', 'pan-human-azimuth'];

// Default values
export const DEFAULT_MAX_CONCURRENCY = 2;
export const DEFAULT_PYTHON_LOG_LEVEL = 40; // Error level
export const DEFAULT_CACHE_DIR = './tmp';
export const DEFAULT_DATASET_HANDLERS = ['hubmap', 'sennet', 'gtex', 'cellxgene', 'ts2'];
export const DEFAULT_DATASET_LIST = 'listing.csv';
export const DEFAULT_DATASET_MIN_CELL_COUNT = 100;
export const DEFAULT_DATASET_DATA_MAX_SIZE = 10 * 1024 * 1024 * 1024; // 10GiB
