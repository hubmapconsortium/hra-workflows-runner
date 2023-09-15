import { opendir } from 'node:fs/promises';
import { join } from 'node:path';
import { env } from 'node:process';
import { findGenerator } from './job-generators/index.js';

const OUTPUT_DIR = env['DATA_REPO_DIR'] || '';

async function getDirectories() {
  try {
    const dir = await opendir(OUTPUT_DIR);
    const dirs = [];
    for await (const entry of dir) {
      if (entry.isDirectory()) {
        if (!entry.path) {
          entry.path = join(dir.path, entry.name);
        }

        dirs.push(entry);
      }
    }

    return dirs;
  } catch (error) {
    console.error(error);
    return [];
  }
}

function groupByGenerator(dirs) {
  const mapping = new Map();
  for (const dir of dirs) {
    const gen = findGenerator(dir);
    const group = mapping.get(gen);
    group ? group.push(dir) : mapping.set(gen, [dir]);
  }

  const items = Array.from(mapping.entries());
  return items.filter(([gen]) => gen !== undefined);
}

async function main() {
  const dirs = await getDirectories();
  for (const [gen, entries] of groupByGenerator(dirs)) {
    await gen.generate(entries);
  }
}

main();
