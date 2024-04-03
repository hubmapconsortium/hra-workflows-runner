/**
 * Base class for custom errors.
 * Ensures that the name property is set correctly.
 */
export class BaseError extends Error {
  constructor(msg) {
    super(msg);
    this.name = this.constructor.name;
  }
}

/**
 * Error indicating an unknown organ code.
 * Stops processing of the associated dataset.
 */
export class UnknownOrganError extends BaseError {
  constructor(code) {
    super(`Unknown organ code '${code}'`);
    this.code = code;
  }
}

/**
 * Formats size as a human readable string
 * Adapted from: https://stackoverflow.com/a/1094933
 *
 * @param {number} size Size in bytes
 */
function formatSize(size) {
  const units = ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi'];
  for (const unit of units) {
    if (size < 1024) {
      return `${size.toFixed(2)}${unit}B`;
    }
    size /= 1024;
  }
  return `${size.toFixed(3)}YiB`;
}

/**
 * Error indicating that the data.h5ad file is too large to be supported.
 * Stops processing of the associated dataset.
 */
export class DataTooLargeError extends BaseError {
  constructor(size) {
    super(`Data size is to large (${formatSize(size)})`);
    this.size = size;
  }
}
