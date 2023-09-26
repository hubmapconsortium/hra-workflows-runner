import { resolve } from 'node:path';
import {
  CACHE_DIR,
  DATASET_HANDLERS,
  DATA_REPO_DIR,
  DEFAULT_MAX_CONCURRENCY,
  FORCE,
  MAX_CONCURRENCY,
  OUTPUT_DIR,
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

    case DATASET_HANDLERS:
      return value.split(/[\s,;]/g);

    case OUTPUT_DIR:
    case DATA_REPO_DIR:
    case CACHE_DIR:
    case SRC_DIR:
      return resolve(value);

    default:
      return value;
  }
}
