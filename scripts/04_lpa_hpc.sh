#!/bin/bash
#$ -cwd
#$ -N lpa_pub
#$ -l h_rt=4:00:00
#$ -l h_data=32G
#$ -j y
#$ -o logs/lpa_pub_$JOB_ID.log
#$ -m bea
#$ -M clarefmccann@g.ucla.edu

. /u/local/Modules/default/init/bash

echo "============================================"
echo "LPA job"
echo "Nodes: $HOSTNAME"
echo "Start: $(date)"
echo "============================================"

module load R/4.2.2
module load gcc/10.2.0

export DATA_DIR="/u/project/silvers/data/ABCD/ABCD-release-6.0/cfm/physical-health/puberty"
export OUT_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"

mkdir -p "${OUT_DIR}/lpa"

export R_MAX_VSIZE=28Gb

SCRIPT_DIR="$SGE_O_WORKDIR"
Rscript "${SCRIPT_DIR}/04_lpa.R"

echo "============================================"
echo "End: $(date)"
echo "============================================"
