#!/bin/bash
source constants.sh
shopt -s extglob
set -ev

CWL_PIPELINE="https://raw.githubusercontent.com/hubmapconsortium/hra-workflows/main/pipeline.cwl"
CWL_OPTS=()

if [[ $RUNNER == "slurm" ]]; then
  CWL_OPTS+=(--singularity)
  SBATCH_FILE="$OUTPUT_DIR/sbatch.sh"

  : >$SBATCH_FILE # Create/truncate the file
  chmod +x $SBATCH_FILE
  echo "#!/bin/bash" >>$SBATCH_FILE
  echo "source constants.sh" >>$SBATCH_FILE
elif [[ $RUNNER == "singularity" ]]; then
  CWL_OPTS+=(--singularity)
fi

if [[ -n $TEMP ]]; then
  CWL_OPTS+=(--tmpdir-prefix $TEMP)
fi

if [[ -d $SIF_CACHE_DIR && ${CWL_OPTS[@]} =~ "singularity" ]]; then
  LINK_SIF_CACHE=true
fi

link_sif_cache() {
  for SIF in $SIF_CACHE_DIR/*.sif; do
    FILE=$(basename $SIF)
    ln -sf $SIF "$1/$FILE"
  done
}

for DIR in $(node $SRC_DIR/list-downloaded-dirs.js); do
  if [[ -n $LINK_SIF_CACHE ]]; then
    link_sif_cache $DIR
  fi

  if [[ $RUNNER == "slurm" ]]; then
    echo "pushd $DIR" >>$SBATCH_FILE
    echo "sbatch $PROJECT_DIR/src/slurm/slurm-annotate.sh" >>$SBATCH_FILE
    echo "popd" >>$SBATCH_FILE
  else
    pushd $DIR
    cwl-runner ${CWL_OPTS[@]} $CWL_PIPELINE job.json
    popd
  fi
done
