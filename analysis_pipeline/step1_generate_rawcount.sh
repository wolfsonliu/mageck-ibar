#! /usr/bin/bash

source step0_preparation.sh

for label in ${LABELS_STEP1[*]}; do
    fq1=${label}.fastq
    fq1gz=${fq1}.gz
    fq2=${label}.fastq
    fq2gz=${fq2}.gz

    ####################

    echo "$(date) unzip ${fq1gz}"
    gunzip ${DIR_STEP0_FASTQ}/${fq1gz}
    echo "$(date) unzip ${fq2gz}"
    gunzip ${DIR_STEP0_FASTQ}/${fq2gz}

    count_sgrna_with_barcode -a "ACCG([ATGC]{20})GTTT[ATGC]{1,35}TGGA([ATCG]{4,6})AACA" \
        -p "ACCG([ATGC]{20})GTTT" \
        -b "TGGA([ATCG]{6})AACA" \
        -f ${DIR_STEP0_FASTQ}/${fq1} -r ${DIR_STEP0_FASTQ}/${fq2} \
        > ${DIR_STEP1_RAWCOUNT}/${label}.rawcount

    echo "$(date) zip ${fq1}"
    gzip ${DIR_STEP0_FASTQ}/${fq1}
    echo "$(date) zip ${fq2}"
    gzip ${DIR_STEP0_FASTQ}/${fq2}

done
####################
# echo "END"
####################
