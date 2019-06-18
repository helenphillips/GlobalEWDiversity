#!/bin/bash
 
#$ -S /bin/bash
#$ -N ESALayer

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -j y

#$ -l h_rt=24:00:00
#$ -l h_vmem=40G,highmem
 
#$ -binding linear:1

set -x

export LANG=en_US.UTF8

output_dir=/work/$USER/$JOB_NAME/$JOB_ID
mkdir -p $output_dir
esa_dir=/data/idiv_sdiv/sworm/esa
soil_dir=/data/idiv_sdiv/sworm/GlobalLayers

module load R
 
Rscript /home/phillips/PrepareESALayer.R $esa_dir $output_dir $soil_dir