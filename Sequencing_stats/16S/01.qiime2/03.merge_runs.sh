#THE STEPS OF THE MERGE:
#1) import all seq runs separately into qiime from manifests
#2) run dada2 on each of these demux.qza files
#3) continuing below, merge all the feature tables and rep seqs on local machine

#merge the otu tables
qiime feature-table merge \
  --i-tables 02_DADA/*table-dada2.qza \
  --o-merged-table artifacts/wd_merged_otu_all.qza

#merge the rep-seqs
qiime feature-table merge-seqs \
  --i-data 02_DADA/*seqs-dada2.qza \
  --o-merged-data artifacts/wd_merged_seqs_all.qza

#export the merged seqs file
qiime tools export \
  --input-path wd_merged_seqs_all.qza \
  --output-path wd_merged_seqs

#export the merged otu table
qiime tools export \
  --input-path wd_merged_otu_all.qza \
  --output-path wd_merged_otu

#make a phylogenetic tree
qiime alignment mafft \
  --i-sequences wd_merged_seqs_all.qza \
  --o-alignment wd_aligned_seqs.qza

qiime phylogeny fasttree \
  --i-alignment wd_aligned_seqs.qza \
  --o-tree wd_unrooted.qza

qiime phylogeny midpoint-root \
  --i-tree wd_unrooted.qza \
  --o-rooted-tree wd_rooted.qza

qiime tools export \
  --input-path wd_rooted.qza \
  --output-path wd_rooted

#convert merged otu table to .biom file
biom convert -i wd_merged_otu.biom -o wd_merged_otu.txt --to-tsv

#to get a taxonomy table, which you will want, open up the merged otu table 
#or rep seqs in the qiime2 viewer and export the taxonomy to tsv
#then modify in excel so you can import into phyloseq in R
