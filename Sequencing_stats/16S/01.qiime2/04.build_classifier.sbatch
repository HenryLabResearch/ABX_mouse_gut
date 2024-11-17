#!/bin/bash
##SBATCH --job-name=buildClassifier
##SBATCH --output=buildClassifier.out
#SBATCH --error=buildClassifier.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load python
conda activate qiime2-2022.2

#step 1: extract reads from s99 database
qiime feature-classifier extract-reads \
  --i-sequences /DIRECTORY/qiime2/S99_otus.qza \
  --p-f-primer AYTGGGYDTAAAGNG \
  --p-r-primer CCGTCAATTYHTTTRAGT \
  --o-reads /DIRECTORY/qiime2/S99_563_926_ref_seqs.qza

#step 2: build classifier
#input reference reads and taxonomy available at https://www.arb-silva.de/download/arb-files/
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /DIRECTORY/qiime2/S99_563_926_ref_seqs.qza \
  --i-reference-taxonomy /DIRECTORY/qiime2/S99ref-taxonomy.qza \
  --p-classify--chunk-size 1000 \
  --o-classifier /DIRECTORY/qiime2/S99_563_926_classifier.qza




