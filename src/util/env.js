import { FORCE } from './constants.js';

export function defaultEnvReviver(key, value) {
  switch (key) {
    case FORCE:
      return value !== '' && value.toLowerCase() !== 'false';

    default:
      return value;
  }
}
