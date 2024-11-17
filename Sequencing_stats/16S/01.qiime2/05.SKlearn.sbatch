#!/bin/bash
##SBATCH --job-name=SKlearn
##SBATCH --output=SKlearn.out
#SBATCH --error=SKlearn.err
#SBATCH --time=10:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=50GB
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.bashrc
module load python
conda activate qiime2-2022.2

qiime feature-classifier classify-sklearn \
  --i-classifier /DIRECTORY/S99_563_926_classifier.qza \
  --i-reads /DIRECTORY/wd_merged_seqs.qza  \
  --p-n-jobs 4 \
  --o-classification /DIRECTORY/wd_merged_s99_taxonomy.qza

qiime metadata tabulate \
  --m-input-file /DIRECTORY/wd_merged_s99_taxonomy.qza \
  --o-visualization /DIRECTORY/wd_merged_s99_taxonomy.qzv

#then download the metadata .tsv file from qiime2 viewer and reformat in excel
#so that it can be imported into phyloseq



