#! /usr/bin/bash

BASE_DIR=/gpfs/share/home/1501111485/Project/Little/20190102
CODE_DIR=${BASE_DIR}/process/p5_plotresult
FIG_DIR=${CODE_DIR}/fig
DATA_DIR=${BASE_DIR}/data
IBAR_DIR=${DATA_DIR}/ibar

for x in A D1 M U1 U4; do
    Rscript ${CODE_DIR}/p5_plotresult.r \
        ${IBAR_DIR}/${x}/iBAR_${x} ${FIG_DIR}/iBAR_${x}
done
