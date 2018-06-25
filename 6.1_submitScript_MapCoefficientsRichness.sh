#!/bin/bash
 
#$ -S /bin/bash
#$ -N MapCoefficients_Richness

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -j y

#$ -l h_rt=100:00:00
#$ -l h_vmem=50G,highmem
 
#$ -binding linear:1

set -x

export LANG=en_US.UTF8

output_dir=/work/$USER/$JOB_NAME/$JOB_ID
mkdir -p $output_dir
GLs_dir=/data/idiv_sdiv/sworm/ProcessGLs
Models_dir=/data/idiv_sdiv/sworm/Models
region="$@"

module load R
 
Rscript /home/phillips/6.1_MapCoefficients_spRichness.R $GLs_dir $Models_dir $output_dir $region