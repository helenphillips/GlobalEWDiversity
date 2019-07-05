#!/bin/bash
 
#$ -S /bin/bash
#$ -N AbundanceBiomassModel

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -j y

#$ -l h_rt=100:00:00
#$ -l h_vmem=100G
 
#$ -binding linear:1

set -x

export LANG=en_US.UTF8

output_dir=/work/$USER/$JOB_NAME/$JOB_ID
mkdir -p $output_dir
data_in=/data/idiv_sdiv/sworm/results/


module load R
 
Rscript /home/phillips/8.2_ModellingAbundanceAndBiomass.R $output_dir $data_in 