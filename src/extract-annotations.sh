#!/bin/bash

if [ -e .venv ]; then
  source .venv/bin/activate
fi
source constants.sh
shopt -s extglob

H5AD=$1
CSV=$(dirname $H5AD)/annotations.csv.gz

python -c "import anndata; ad=anndata.read_h5ad('${H5AD}', backed='r'); ad.obs.to_csv('${CSV}', compression='gzip'); print('finished')"
