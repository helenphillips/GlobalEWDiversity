#!/bin/bash
 
#$ -S /bin/bash
#$ -N GlobalDataLayers

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -j y

#$ -l h_rt=100:00:00
#$ -l h_vmem=50G,highmem
 
#$ -binding linear:1

set -x

export LANG=en_US.UTF8

output_dir=/work/$USER/$JOB_NAME/$JOB_ID
mkdir -p $output_dir
data_dir=/data/idiv_sdiv/sworm/GlobalLayers
date="2019-11-25"
site_dir=/data/idiv_sdiv/sworm
module load R
 
Rscript /home/phillips/PrepareGlobalLayers_abundance.R $data_dir $output_dir $date $site_dir