#!/bin/bash
##SBATCH --job-name=samtools_samtobam_CO
##SBATCH --output=samtools.out
#SBATCH --error=samtools.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=2200
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.bashrc
module load python
conda activate samtools


for sample in `cat LIST`;
do
  	samtools view -F 4 -bS -@ 20 05_MAPPING/$sample.sam > 05_MAPPING/$sample.bam
        samtools sort -@ 20 05_MAPPING/$sample.bam -o 05_MAPPING/$sample.sorted.bam
        samtools index -@ 20 05_MAPPING/$sample.sorted.bam
done
