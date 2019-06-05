#! /usr/bin/bash
#SBATCH --job-name=B_p_ibartest
#SBATCH --output=/gpfs/share/home/1501111485/log/%x.%j_%A_%a.%N.out
#SBATCH --error=/gpfs/share/home/1501111485/log/%x.%j_%A_%a.%N.err
#SBATCH --partition=C032M0128G
#SBATCH --qos=low
#SBATCH --get-user-env
#SBATCH -n 2
#SBATCH --cpu-freq=high
#SBATCH -A hpc1706879002
#SBATCH --mail-type=end
#SBATCH --mail-user=1501111485@pku.edu.cn
#SBATCH --time=120:00:00

module load anaconda/3-4.4.0.1
# source activate ibar
source activate ibar

labels=(A D1 M U1 U4)

label=${labels[$SLURM_ARRAY_TASK_ID]}

####################

BASE_DIR=/gpfs/share/home/1501111485/Project/Little/20190102
CODE_DIR=${BASE_DIR}/process/p2_librarycount
DATA_DIR=${BASE_DIR}/data
FQ_DIR=${DATA_DIR}/raw
RAWCOUNT_DIR=${DATA_DIR}/rawcount
COUNT_DIR=${DATA_DIR}/count
IBAR_DIR=${DATA_DIR}/ibar

####################

echo "$(date) Start: "${label}

if [ ! -e ${IBAR_DIR}/${label} ]; then
    mkdir -p ${IBAR_DIR}/${label}
fi

cp ${COUNT_DIR}/${label}.count.csv ${IBAR_DIR}/${label}/

cd ${IBAR_DIR}/${label}

mageck-ibar -i ${IBAR_DIR}/${label}/${label}.count.csv \
    -b --col-gene gene --col-guide guide --col-barcode barcode \
    -c Ctr_R1 Ctr_R2 -t ${label}_R1 ${label}_R2 \
    --largerthan 10 \
    -o ${IBAR_DIR}/${label}/iBAR_${label}

echo "$(date) FINISH"

####################
source deactivate
module unload anaconda/3-4.4.0.1
