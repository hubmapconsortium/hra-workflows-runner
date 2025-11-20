import { execFile as callbackExecFile } from 'node:child_process';
import { promisify } from 'node:util';
import { Config } from '../util/config.js';
import { getSrcFilePath } from '../util/paths.js';

const execFile = promisify(callbackExecFile);

/**
 * Infers RNA source (cell vs nucleus) from an h5ad file using the simplified
 * spliced/unspliced heuristic implemented in the Python helper.
 *
 * @param {string} h5adPath Path to the h5ad file
 * @param {Config} config Configuration object
 * @param {number} [unsplicedThreshold=0.30] Threshold for unspliced fraction to call nucleus
 * @returns {Promise<object>} Inference results with verdict, evidence, and metrics
 */
export async function InferRnaSourceFromH5ad(h5adPath, config, unsplicedThreshold = 0.30) {

  const scriptPath = getSrcFilePath(config, 'xconsortia', 'infer_prep_from_h5ad.py');

  try {
    const args = [
      scriptPath,
      h5adPath,
      '--unspliced-threshold',
      unsplicedThreshold.toString(),
    ];
    const { stdout } = await execFile('python3', args);

    // Parse the JSON output from stdout
    const result = JSON.parse(stdout);
    return result;
  } catch (error) {

    console.warn(`Failed to infer prep from h5ad file ${h5adPath}: ${error.message}`);
    return {
      input_file: h5adPath,
      verdict: 'error',
      error: error.message,
      evidence: [`Failed to run inference: ${error.message}`],
      metrics: {},
    };
  }
}

