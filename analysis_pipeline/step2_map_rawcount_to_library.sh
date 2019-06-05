#! /usr/bin/bash

source step0_preparation.sh

for label in ${LABELS_STEP2[*]}; do
    python3 ${DIR_STEP2}/makeinput.py --reference ${FILE_STEP0_LIBRARY} \
        --controllabel ${LABEL_STEP2_CTRL}_1  ${LABEL_STEP2_CTRL}_2 \
        --controlinput ${DIR_STEP2_RAWCOUNT}/${LABEL_STEP2_CTRL}_1.rawcount \
        ${DIR_STEP2_RAWCOUNT}/${LABEL_STEP2_CTRL}_2.rawcount \
        --treatlabel ${label}_1 ${label}_2 \
        --treatinput ${DIR_STEP2_RAWCOUNT}/${label}_1.rawcount \
        ${DIR_STEP2_RAWCOUNT}/${label}_2.rawcount \
        --output ${DIR_STEP2_COUNT}/${label}.count.csv
done
####################
# echo "END"
####################
