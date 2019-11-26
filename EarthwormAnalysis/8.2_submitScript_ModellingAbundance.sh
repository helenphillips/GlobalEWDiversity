#!/bin/bash
 
#$ -S /bin/bash
#$ -N AbundanceModel

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -j y

#$ -l h_rt=100:00:00
#$ -l h_vmem=10G
 
#$ -binding linear:1

set -x

export LANG=en_US.UTF8

output_dir=/work/$USER/$JOB_NAME/$JOB_ID
mkdir -p $output_dir
data_in=/data/idiv_sdiv/sworm
date="2019-11-25"
functions=/work/$USER/Functions


module load R
 
Rscript /home/phillips/8.2_ModellingAbundance.R $output_dir $data_in $date $functions