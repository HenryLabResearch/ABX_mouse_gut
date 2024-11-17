#!/bin/bash
##SBATCH --job-name=bowtie2-build
##SBATCH --output=bowtie-build.out
#SBATCH --error=bowtie.err
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

        bowtie2-build 04_CONTIGS/contigs.fixed.fa 05_MAPPING/contigs --threads 22 -v

