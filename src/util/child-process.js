import { spawn } from 'node:child_process';
import { createWriteStream } from 'node:fs';
import { finished } from 'node:stream/promises';

const DEFAULT_STDOUT_MAX_BYTES = 1024 * 1024;
const DEFAULT_STDERR_MAX_BYTES = 64 * 1024;

/**
 * Runs a child process with bounded output capture to avoid large string growth.
 *
 * @param {string} command Command to execute
 * @param {string[]} args Command arguments
 * @param {{
 *   captureStdout?: boolean,
 *   captureStdoutLimit?: number,
 *   captureStderrLimit?: number,
 *   stdoutFile?: string,
 *   stderrFile?: string,
 *   onStdoutChunk?: (chunk: Buffer) => void,
 *   onStderrChunk?: (chunk: Buffer) => void,
 *   spawnOptions?: import('node:child_process').SpawnOptions,
 * }} [options]
 * @returns {Promise<{stdout: string, stderr: string}>}
 */
export async function spawnProcess(command, args, options = {}) {
  const {
    captureStdout = true,
    captureStdoutLimit = DEFAULT_STDOUT_MAX_BYTES,
    captureStderrLimit = DEFAULT_STDERR_MAX_BYTES,
    stdoutFile,
    stderrFile,
    onStdoutChunk,
    onStderrChunk,
    spawnOptions = {},
  } = options;

  return new Promise((resolve, reject) => {
    const child = spawn(command, args, {
      stdio: ['ignore', 'pipe', 'pipe'],
      ...spawnOptions,
    });
    const stdoutStream = stdoutFile ? createWriteStream(stdoutFile) : undefined;
    const stderrStream = stderrFile ? createWriteStream(stderrFile) : undefined;
    let stdout = '';
    let stderr = '';
    let stdoutBytes = 0;
    let stderrBytes = 0;
    let settled = false;
    /** @type {Error | undefined} */
    let processError;

    const handleProcessError = (error) => {
      if (processError) {
        return;
      }

      processError = error;
      child.kill();
    };

    const settle = (finalizer) => {
      if (settled) {
        return;
      }

      settled = true;
      stdoutStream?.end();
      stderrStream?.end();

      Promise.allSettled(
        [stdoutStream, stderrStream]
          .filter((stream) => stream !== undefined)
          .map((stream) => finished(stream))
      ).finally(finalizer);
    };

    stdoutStream?.on('error', handleProcessError);
    stderrStream?.on('error', handleProcessError);

    const appendChunk = (text, chunk, currentBytes, limit, streamName) => {
      if (limit <= 0) {
        return { text, bytes: currentBytes };
      }

      const chunkBytes = chunk.byteLength;
      const nextBytes = currentBytes + chunkBytes;
      if (nextBytes > limit) {
        handleProcessError(new Error(`${command} ${streamName} exceeded ${limit} bytes`));
        return { text, bytes: currentBytes };
      }

      return {
        text: text + chunk.toString(),
        bytes: nextBytes,
      };
    };

    child.stdout.on('data', (chunk) => {
      stdoutStream?.write(chunk);
      onStdoutChunk?.(chunk);
      if (!captureStdout || processError) {
        return;
      }

      const next = appendChunk(stdout, chunk, stdoutBytes, captureStdoutLimit, 'stdout');
      stdout = next.text;
      stdoutBytes = next.bytes;
    });

    child.stderr.on('data', (chunk) => {
      stderrStream?.write(chunk);
      onStderrChunk?.(chunk);
      if (processError) {
        return;
      }

      const next = appendChunk(stderr, chunk, stderrBytes, captureStderrLimit, 'stderr');
      stderr = next.text;
      stderrBytes = next.bytes;
    });

    child.on('error', (error) => {
      settle(() => reject(error));
    });
    child.on('close', (code, signal) => {
      if (processError) {
        settle(() => reject(processError));
        return;
      }

      if (code === 0) {
        settle(() => resolve({ stdout, stderr }));
        return;
      }

      const reason = signal ? `signal ${signal}` : `exit code ${code}`;
      const error = new Error(`${command} failed with ${reason}${stderr ? `: ${stderr.trim()}` : ''}`);
      settle(() => reject(error));
    });
  });
}