#!/bin/bash
##SBATCH --job-name=MERGE_PROFILES
##SBATCH --output=mergeProfiles.out
#SBATCH --error=mergeProfiles.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=22
#SBATCH --mem-per-cpu=2100
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
 
source ~/.bashrc
module load python
conda activate anvio-7
 
        anvi-merge 05_MAPPING/*/PROFILE.db -c 04_CONTIGS/contigs.db -o 06_MERGED_PROFILES

