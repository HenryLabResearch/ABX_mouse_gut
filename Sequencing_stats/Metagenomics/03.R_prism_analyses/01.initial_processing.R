library(tidyverse)

setwd("~/Documents/Chang-Bergelson Lab/HF Diet Abx/Aim 1/Metagenomics/comm_analysis/")

#this is the gene call coverages output file from anvio
cov = read.csv("wd_all_cov.csv")

#this file maps the gene calls to KEGG Orthologs (KOs)
ko_genes = read_table("wd_kofams.txt")
ko_genes = dplyr::select(ko_genes, c("gene_callers_id","accession")) %>%
  rename(accession = "KO")

#filter out gene calls that don't map to KEGG
good = cov$gene_callers_id[cov$gene_callers_id %in% ko_genes$gene_callers_id]
cov1 = cov[cov$gene_callers_id %in% good,]

#map the gene call coverages to their KO
cov2 = left_join(cov1, ko_genes, multiple="all")
cov2 = cov2[,-1] #get rid of gene calls column

#add up gene calls that map to the same KO
ko_sum = plyr::ddply(cov2, "KO", plyr::numcolwise(sum))
row.names(ko_sum) = ko_sum$KO

#write.csv(ko_sum, "counts_ko_sum.csv")

