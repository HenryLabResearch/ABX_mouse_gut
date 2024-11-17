library(tidyverse)
library(phyloseq)
library(microbiome)
library(picante)

setwd("/DIRECTORY/")

meta = read.csv("meta.csv")
row.names(meta) = meta$SeqID


### KO-level first ###

ko <-read.csv("counts_ko_sum.csv",row.names=1) #this file has KO-level abundances
ko1 = ko[colnames(ko) %in% meta$SeqID]

ko.m = as.matrix(ko1)

OTU = otu_table(ko.m, taxa_are_rows = TRUE)
META = sample_data(meta)
ps = phyloseq(OTU,META)

ko.r = microbiome::richness(ps)
ko.r$SeqID = rownames(ko.r)

ko.r.meta = left_join(ko.r, meta)
ko.r.meta$Treatment = factor(ko.r.meta$Treatment, levels=c("WD-ABX","RC-ABX"))

#gather the pre-abx timepoint so we can normalize to that
ko.r.bl = subset(ko.r.meta, Rec_day_adj=="-3")
ko.r.bl = dplyr::select(ko.r.bl, c(observed, Mouse))

#pre_rich is the richness at the pre-abx timepoint
colnames(ko.r.bl)[colnames(ko.r.bl)=="observed"] = "pre_rich"
ko.r.bl$pre = "pre"

ko.r.meta2 = left_join(ko.r.meta,ko.r.bl)

#richness percentage = observed at d14 divided by d-3
ko.r.meta2$rich_perc = ko.r.meta2$observed/ko.r.meta2$pre_rich

ko.rich.perc.summ = ko.r.meta2 %>% group_by(Rec_day_adj, Treatment) %>% 
  summarise(
    mean_richness = mean(rich_perc),
    sd = sd(rich_perc, na.rm=TRUE))

ko.rich.perc.summ$level = "KO"


### Gene Call level ###

genes <-read.csv("wd_all_cov.csv",row.names=1)
genes1 = genes[colnames(genes) %in% meta$SeqID]

genes.m = as.matrix(genes1)

OTU = otu_table(genes.m, taxa_are_rows = TRUE)
META = sample_data(meta)
ps = phyloseq(OTU,META)

gene.r = microbiome::richness(ps)
gene.r$SeqID= rownames(gene.r)

genes.meta = left_join(gene.r, meta)
genes.meta$Treatment = factor(genes.meta$Treatment, levels=c("WD-ABX","RC-ABX"))

#pre-abx values for comparison
genes.bl = subset(genes.meta, Rec_day_adj=="-3")
genes.bl = dplyr::select(genes.bl, c(observed, Mouse))
colnames(genes.bl)[colnames(genes.bl)=="observed"] = "pre_rich"
genes.bl$pre = "pre"

genes.meta2 = left_join(genes.meta,genes.bl)
genes.meta2$rich_perc = genes.meta2$observed/genes.meta2$pre_rich

genes.rich.perc.summ = genes.meta2 %>% group_by(Rec_day_adj, Treatment) %>% 
  summarise(
    mean_richness = mean(rich_perc),
    sd = sd(rich_perc, na.rm=TRUE))

genes.rich.perc.summ$level = "Gene Call"

### KEGG Category level

ko1$KO = rownames(ko1)

ko.levels = read.csv("ko.all.csv")
kcats = dplyr::select(ko.levels, c(KEGG.category, KO))

kcat = left_join(ko1, kcats, multiple="all")
kcat2 = select(kcat, -KO)
kcat3 = kcat2[!is.na(kcat2$KEGG.category),] #drop the NA

#sum up all values for same KCat
kcat_sum = plyr::ddply(kcat3, "KEGG.category", plyr::numcolwise(sum))

rownames(kcat_sum) = kcat_sum$KEGG.category
kcat_sum = kcat_sum[,-1]

kcat.m = as.matrix(kcat_sum)

OTU = otu_table(kcat.m, taxa_are_rows = TRUE)
META = sample_data(meta)
ps = phyloseq(OTU,META)

kc.r = microbiome::richness(ps)
kc.r$SeqID= rownames(kc.r)

kc.r.meta = left_join(kc.r, meta)

kc.r.meta$Treatment = factor(kc.r.meta$Treatment, levels=c("WD-ABX","RC-ABX"))


kc.r.bl = subset(kc.r.meta, Rec_day_adj=="-3")
kc.r.bl = dplyr::select(kc.r.bl, c(observed, Mouse))
colnames(kc.r.bl)[colnames(kc.r.bl)=="observed"] = "pre_rich"
kc.r.bl$pre = "pre"

kc.r.meta2 = left_join(kc.r.meta,kc.r.bl)
kc.r.meta2$rich_perc = kc.r.meta2$observed/kc.r.meta2$pre_rich

kcat.rich.perc.summ = kc.r.meta2 %>% group_by(Rec_day_adj, Treatment) %>% 
  summarise(
    mean_richness = mean(rich_perc),
    sd = sd(rich_perc, na.rm=TRUE))

kcat.rich.perc.summ$level = "KEGG Category"


######### put the summaries together
richness_all.levels = rbind(ko.rich.perc.summ, genes.rich.perc.summ, kcat.rich.perc.summ)
#write.csv(richness_all.levels, "richness_all_levels_percent.csv",row.names=FALSE)

# ^^^ this was plotted and further analyzed in 04.func_richness_redund.pzfx
