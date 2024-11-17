#!/bin/bash
##SBATCH --job-name=megahit.CO
##SBATCH --output=megahit.out
#SBATCH --error=megahit.err
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=22
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.bashrc
module load python
conda activate anvio-7

R1s=`ls 02_QC/*QUALITY_PASSED_R1* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`

R2s=`ls 02_QC/*QUALITY_PASSED_R2* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`

megahit -1 $R1s -2 $R2s --min-contig-len 2000 -m 0.85 -o 03_COASSEMBLY/ -t 22