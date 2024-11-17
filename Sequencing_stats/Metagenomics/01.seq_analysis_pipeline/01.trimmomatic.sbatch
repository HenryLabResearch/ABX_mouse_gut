#!/bin/bash
##SBATCH --job-name=TRIM_QC1
##SBATCH --output=TRIM.out
#SBATCH --error=TRIM.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mkennedy21@uchicago.edu

source ~/.bashrc
module load java
module load python
conda activate TRIMMOMATIC

for sample in `cat LIST3`;
do
        java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -phred33 00_RAW/$sample.R1.fastq.gz 00_RAW/$sample.R2.fastq.gz 01_TRIMMED/$sample.trimmed.R1.fastq.gz 01_TRIMMED/$sample.unpaired.R1.fastq.gz 01_TRIMMED/$sample.trimmed.R2.fastq.gz 01_TRIMMED/$sample.unpaired.R2.fastq.gz ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
