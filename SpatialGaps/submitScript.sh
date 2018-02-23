#!/bin/bash
 
#$ -S /bin/bash
#$ -N accessibility

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -j y

#$ -l h_rt=8:00:00
#$ -l h_vmem=100G,highmem
 
#$ -binding linear:1

set -x

export LANG=en_US.UTF8

output_dir=/work/$USER/$JOB_NAME/$JOB_ID
mkdir -p $output_dir
data_dir=/data/idiv_sdiv

module load R
 
Rscript /home/phillips/AccessibilityAndClimate.R $data_dir $output_dir