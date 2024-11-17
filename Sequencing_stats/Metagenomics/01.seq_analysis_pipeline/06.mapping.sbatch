#!/bin/bash
##SBATCH --job-name=bowtie2-mapping_CO
##SBATCH --output=mapping.out
#SBATCH --error=mapping.err
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

        bowtie2 --threads 22 -x 05_MAPPING/contigs \
        -1 02_QC/$sample-QUALITY_PASSED_R1.fastq -2 02_QC/$sample-QUALITY_PASSED_R2.fastq \
        -S 05_MAPPING/$sample.sam

done
