#! /usr/bin/bash

DIR_BASE=./analysis

####################
# Step0: Prepare for analysis
DIR_STEP0=${DIR_BASE}/step0_preparation
DIR_STEP0_FASTQ=${DIR_STEP0}/fastq
FILE_STEP0_LIBRARY=${DIR_STEP0}/library.txt

####################
# Step1: Generate rawcounts from fastq files
DIR_STEP1=${DIR_BASE}/step1_generate_rawcount
DIR_STEP1_RAWCOUNT=${DIR_STEP1}/rawcount
LABELS_STEP1=(Ctrl_1 Ctrl_2 Exp_1 Exp_2)

####################
# Step2: Merge rawcounts to library file
DIR_STEP2=${DIR_BASE}/step2_map_rawcount_to_library
DIR_STEP2_COUNT=${DIR_STEP2}/count
LABELS_STEP2=(Exp)
LABEL_STEP2_CTRl=Ctrl

####################
# Step3: Quality check of counts
DIR_STEP3=${DIR_BASE}/step3_quality_check
DIR_STEP3_FIG=${DIR_STEP3}/fig

####################
# Step4: Run mageck-ibar analysis
DIR_STEP4=${DIR_BASE}/step4_mageck-ibar_analysis
DIR_STEP4_IBAR=${DIR_STEP4}/ibar

####################
# Step5: Analysis result visualization
DIR_STEP5=${DIR_BASE}/step5_analysis_result
DIR_STEP5_FIG=${DIR_STEP5}/fig
####################
# echo "END"
####################
