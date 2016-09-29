#!/bin/sh
#PBS -N parallel_combine_target_databases
#PBS -l nodes=1:ppn=8,walltime=240:00:00

cd $PBS_O_WORKDIR   
module load r/3.1.1
R CMD BATCH parallel_combine_databases.R
