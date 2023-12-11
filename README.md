# hra-workflows-runner
Code to run hra-workflows locally or using Slurm.

## Requirements
- [cwl-runner](https://github.com/common-workflow-language/cwltool)
- Docker (only when running locally) **OR**
- Singularity (sometimes called Apptainer)

## Running the pipeline

### 1. Configure the environment
Copy `env.example.sh` into `env.sh` and then update the configuration with your own settings.

### 2. Run scripts
If using singularity rather than docker it is highly recommended that you configure the `SIF_CACHE_DIR` option and run `./scripts/00x-build-containers.sh` to prebuild the sif containers for faster runs.

Another script that only needs to be run once is `05x-download-models.sh`. This script will download the model files that were to large to embed directly in the algorithm's container file.

The next step depends on whether you are running the code locally or in an environment using Slurm.

#### Locally
Run `run.sh`. This script will run the entire pipeline from start to finish. Note that the annotation steps will run sequentially for all datasets and may therefore take a lot of time when processing large amounts of datasets.

#### Slurm
Running on Slurm requires some manual work due to difficulties of running nested singularity containers. Start by running `x-run.sh`. This script will schedule a job that will run every step up to but not including annotating the datasets. Once the `x-run.sh` job has finished the annotation step can be started using `./scripts/30x-annotate.sh`. After all annotations have finished the `40+` scripts have to be run manually. `./scripts/01x-start-container.sh` is provided as a utility to enter an environment where `40+` scripts can be run with all dependencies satisfied.

## Adding new dataset handlers
Adding a new dataset handler requires two steps:

### 1. Implementing interfaces
A dataset handler must implement and export the [DatasetHandler](./src/util/handler.js) interface from it's index.js file. The interface includes a listing generator, downloader, and job generator. It can be useful to refer to one of the existing implementations such as [CellXGene](./src/cellxgene/index.js) or [GTEx](./src/gtex/index.js) when creating a new handler.

#### Note on source location:
The `index.js` file must be placed in a `src/handler-name/` folder to allow the handler to be properly loaded.

### 2. Registering the new handler
The new handler must then be registered by adding it's name in [DEFAULT_DATASET_HANDLERS](./src/util/constants.js). Alternatively handlers can be specified using the `DATASET_HANDLERS` environment variable.
