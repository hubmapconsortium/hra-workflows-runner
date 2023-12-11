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
