import {
  DEFAULT_MAX_CONCURRENCY,
  FORCE,
  MAX_CONCURRENCY,
} from './constants.js';

export function defaultEnvReviver(key, value) {
  switch (key) {
    case FORCE:
      return value !== '' && value.toLowerCase() !== 'false';

    case MAX_CONCURRENCY: {
      const num = Number.parseInt(value);
      return num > 0 ? num : DEFAULT_MAX_CONCURRENCY;
    }

    default:
      return value;
  }
}
