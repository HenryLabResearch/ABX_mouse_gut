#!/bin/bash
##SBATCH --job-name=MERGE
##SBATCH --output=ANVI-MERGE.out
#SBATCH --error=ANVI-MERGE.err
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=500000
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.bashrc
module load python
conda activate anvio-7
  
anvi-export-gene-coverage-and-detection -c contigs.db \
                                        -p 06_MERGED_PROFILES/PROFILE.db \
                                        -O wd_all_cov