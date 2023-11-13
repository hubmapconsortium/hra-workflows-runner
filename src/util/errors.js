export class BaseError extends Error {
  constructor(msg) {
    super(msg);
    this.name = this.constructor.name;
  }
}

export class UnknownOrganError extends BaseError {
  constructor(code) {
    super(`Unknown organ code '${code}'`);
    this.code = code;
  }
}
