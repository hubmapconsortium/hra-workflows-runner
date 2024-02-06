import { readFileSync, writeFileSync } from 'fs';
import { globSync } from 'glob';
import Papa from 'papaparse';
import { basename } from 'path';

const crosswalks = globSync('./crosswalking-tables/*.csv');
const OUTPUT = 'crosswalking-tables/crosswalks.jsonld';

const results = [];
for (const crosswalkFile of crosswalks) {
  const tool = basename(crosswalkFile, '.csv');
  const csvStringRows = readFileSync(crosswalkFile).toString().split('\n');
  const headerRow = csvStringRows.findIndex((s) => s.startsWith('Organ_Level'));
  const csvString = csvStringRows.slice(headerRow).join('\n');
  const { data } = Papa.parse(csvString, { skipEmptyLines: true, header: true });

  for (const row of data) {
    results.push({
      '@id': row.Annotation_Label_ID,
      '@type': 'AnnotationItem',
      label: row.Annotation_Label,
      tool,
      organ_id: row.Organ_ID,
      organ_level: row.Organ_Level,
      cell_id: row.CL_ID,
      cell_label: row.CL_Label,
      cell_match: row.CL_Match,
    });
  }
}

// Write out the new crosswalks.jsonld file
const jsonld = {
  ...JSON.parse(readFileSync('src/summary-context.jsonld')),
  '@graph': results,
};
writeFileSync(OUTPUT, JSON.stringify(jsonld, null, 2));
