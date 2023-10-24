import { resolve } from 'node:path';
import {
  CACHE_DIR,
  DATASETS_DIR,
  DATASET_HANDLERS,
  DATA_REPO_DIR,
  DEFAULT_MAX_CONCURRENCY,
  DEFAULT_PYTHON_LOG_LEVEL,
  FORCE,
  MAX_CONCURRENCY,
  MODELS_DIR,
  OUTPUT_DIR,
  PYTHON_LOG_LEVEL,
  SRC_DIR,
} from './constants.js';

export function defaultEnvReviver(key, value) {
  switch (key) {
    case FORCE:
      return value !== '' && value.toLowerCase() !== 'false';

    case MAX_CONCURRENCY: {
      const num = Number.parseInt(value);
      return num > 0 ? num : DEFAULT_MAX_CONCURRENCY;
    }

    case PYTHON_LOG_LEVEL: {
      const num = Number.parseInt(value);
      return num > 0 ? num : DEFAULT_PYTHON_LOG_LEVEL;
    }

    case DATASET_HANDLERS:
      return value.split(/[\s,;]/g);

    case DATASETS_DIR:
    case OUTPUT_DIR:
    case DATA_REPO_DIR:
    case CACHE_DIR:
    case MODELS_DIR:
    case SRC_DIR:
      return resolve(value);

    default:
      return value;
  }
}
