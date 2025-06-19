#!/bin/bash
#SBATCH --mem=20000
#SBATCH --cpus-per-task=20

echo "My job ran on : $(hostname)"
echo "My job started at : $(date)" 

module load FastQC

tcsh run_fastqc_1.sh 


