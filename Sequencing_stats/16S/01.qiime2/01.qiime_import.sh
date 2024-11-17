#this goes through each manifest and imports them into a qiime .qza object

for file in $(find manifests -name 'manifest*.txt'); do

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $file \
  --output-path ${file}.qza \
  --input-format PairedEndFastqManifestPhred33V2

done
