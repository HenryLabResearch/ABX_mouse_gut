#!/bin/bash
##SBATCH --job-name=DADA2
##SBATCH --output=DADA2.out
#SBATCH --error=DADA2.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.bashrc
module load python
conda activate qiime2-2022.2

for file in $(find artifacts -name '*.qza'); do

filename=$(basename "$file")
basename=$(echo "${filename%.*}")

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $file \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-representative-sequences 02_DADA/${basename}.rep-seqs-dada2.qza \
  --o-table 02_DADA/${basename}.table-dada2.qza \
  --o-denoising-stats 02_DADA/${basename}.stats-dada2.qza \
  --p-n-threads 6 \
  --verbose

#no truncation since median quality score was >30 for all datasets at full length

qiime feature-table summarize \
  --i-table 02_DADA/${basename}.table-dada2.qza \
  --o-visualization 02_DADA/${basename}.table-dada2.qzv

qiime feature-table tabulate-seqs \
  --i-data 02_DADA/${basename}.rep-seqs-dada2.qza \
  --o-visualization 02_DADA/${basename}.rep-seqs-dada2.qzv

qiime metadata tabulate \
  --m-input-file 02_DADA/${basename}.stats-dada2.qza \
  --o-visualization 02_DADA/${basename}.stats.qzv

done
