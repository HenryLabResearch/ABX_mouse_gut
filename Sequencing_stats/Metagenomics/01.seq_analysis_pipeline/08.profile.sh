#!/bin/bash
##SBATCH --job-name=PROFILES_CO
##SBATCH --output=profiles.out
#SBATCH --error=profiles.err
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
 
for sample in `cat LIST`;

do 

        anvi-profile -i 05_MAPPING/$sample.sorted.bam -c 04_CONTIGS/contigs.db -T 22

done