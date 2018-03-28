#!/bin/bash
 
#$ -S /bin/bash
#$ -N GlobalDataLayers

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -j y

#$ -l h_rt=24:00:00
#$ -l h_vmem=50G,highmem
 
#$ -binding linear:1

set -x

export LANG=en_US.UTF8

output_dir=/work/$USER/$JOB_NAME/$JOB_ID
mkdir -p $output_dir
data_dir=/data/idiv_sdiv
date="2018-03-27"
module load R
 
Rscript /home/phillips/PrepareGlobalLayers.R $data_dir $output_dir $date