#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows-runner:main
    dockerOutputDirectory: /output
  NetworkAccess:
    networkAccess: true
  EnvVarRequirement:
    envDef:
      DATASET: $(inputs.dataset)
      VERSION: $(inputs.version)
      HUBMAP_TOKEN: $(inputs.hubmap_token)
      DATASETS_DIR: /output/datasets
      RAW_DATA_DIR: /output/raw-data
  InitialWorkDirRequirement:
    listing:
    - entryname: datasets
      writable: false
      entry: $(inputs.datasets_dir)
    - entryname: raw-data
      writable: true
      entry: $(inputs.rawdata_dir)

stdout: output.txt

inputs:
  datasets_dir:
    type: Directory
  rawdata_dir:
    type: Directory
  dataset:
    type: string
    default: sample
  version:
    type: string
    default: v1
  hubmap_token:
    type: string
    default: ''
  command:
    type: string
    default: ./run.sh

outputs:
  rawdata_out:
    type: Directory
    outputBinding:
      glob: raw-data
  command_out:
    type: stdout

baseCommand: bash
arguments:
  - -c
  - cd /workspace && $(inputs.command)
