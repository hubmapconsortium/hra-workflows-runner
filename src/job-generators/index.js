import * as gtex from './gtex.js';
import * as hubmap from './hubmap.js';

const GENERATORS = [gtex, hubmap];

export function findGenerator(dir) {
  return GENERATORS.find((gen) => gen.supports(dir));
}
