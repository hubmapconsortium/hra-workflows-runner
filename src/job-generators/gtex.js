import { execFile as rawExecFile } from 'node:child_process';
import { writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { promisify } from 'node:util';

const PYTHON_SCRIPT = 'src/job-generators/gtex-get-organ.py';

const execFile = promisify(rawExecFile);

function getJobYaml(organ) {
  return `organ: ${organ}
matrix:
  class: File
  path: data.h5ad

preprocessing:
  geneColumn: counts

algorithms:
#  - azimuth: {}
#    extract:
#      annotationMethod: azimuth
#      cellLabelColumn: azimuth_label
#      geneLabelColumn: gene

  - celltypist: {}
    extract:
      annotationMethod: celltypist
      cellLabelColumn: predicted_labels
      geneLabelColumn: gene

#  - popv:
#      referenceData:
#        class: File
#        path: TS_Lung_filtered.h5ad
#    extract:
#      annotationMethod: popv
#      cellLabelColumn: popv_prediction
#      geneLabelColumn: gene
`;
}

async function getOrgan(dir) {
  try {
    const { stdout } = await execFile('python3', [
      PYTHON_SCRIPT,
      join(dir.path, 'data.h5ad'),
    ]);
    const match = /organ:\s*(.+)/i.exec(stdout);
    return match?.[1];
  } catch (error) {
    console.error(error);
    return undefined;
  }
}

async function generateJob(dir, organ) {
  if (!organ) {
    return;
  }

  try {
    const path = join(dir.path, 'job.yaml');
    await writeFile(path, getJobYaml(organ));
  } catch (error) {
    console.error(error);
  }
}

export function supports(dir) {
  return /^gtex/i.test(dir.name);
}

export async function generate(dirs) {
  for (const dir of dirs) {
    const organ = await getOrgan(dir);
    await generateJob(dir, organ);
  }
}
