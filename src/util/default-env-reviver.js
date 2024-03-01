import { resolve } from 'node:path';
import {
  CACHE_DIR,
  CROSSWALKING_TABLES_DIR,
  DATASETS_DIR,
  DATASET_HANDLERS,
  DATASET_MIN_CELL_COUNT,
  DATA_REPO_DIR,
  DEFAULT_DATASET_MIN_CELL_COUNT,
  DEFAULT_MAX_CONCURRENCY,
  DEFAULT_PYTHON_LOG_LEVEL,
  FORCE,
  MAX_CONCURRENCY,
  MODELS_DIR,
  OUTPUT_DIR,
  PYTHON_LOG_LEVEL,
  SRC_DIR,
} from './constants.js';

/**
 * Default environment reviver used in Config#loadEnv
 *
 * @param {string} key Environment variable name
 * @param {string} value Environment variable value
 * @returns {any} The parsed value
 */
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

    case DATASET_MIN_CELL_COUNT: {
      const num = Number.parseInt(value);
      return num >= 0 ? num : DEFAULT_DATASET_MIN_CELL_COUNT;
    }

    case DATASET_HANDLERS:
      return value.split(/[\s,;]/g);

    case DATASETS_DIR:
    case OUTPUT_DIR:
    case DATA_REPO_DIR:
    case CACHE_DIR:
    case MODELS_DIR:
    case CROSSWALKING_TABLES_DIR:
    case SRC_DIR:
      return resolve(value);

    default:
      return value;
  }
}
