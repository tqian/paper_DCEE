#!/bin/bash

#SBATCH --job-name=distal  ## job name
#SBATCH -p standard               ## account to charge
#SBATCH -A tiancq1_lab              ## use lab partition
#SBATCH --nodes=1           ## use 1 node, don't ask for multiple
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=1    ## number of cores the job needs
#SBATCH --mem-per-cpu=4G    ## ask for 1Gb memory per CPU
#SBATCH --error=output/%x.%A.err    ## Slurm error  file, %x - job name, %A job id
#SBATCH --out=output/%x.%A.out      ## Slurm output file, %x - job name, %A job id

module load R/4.3.3 # load R
Rscript simu-distal.R $1
