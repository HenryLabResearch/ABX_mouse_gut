library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(sva)

setwd("/DIRECTORY/")

ko_sum = read.csv("counts_ko_sum.csv", row.names=1)

meta = read.csv("meta.csv")

meta = subset(meta, Rec_day_adj != 8 & perc_mapping > 0.5)
meta$Rec_day_adj = factor(meta$Rec_day_adj, levels=c("-13","-3","2","4","14","28"))


#subset by diet
rc = subset(meta, Treatment=="RC-ABX")
wd = subset(meta, Treatment=="WD-ABX")


#filter to the samples in the metadata
rc.counts = ko_sum[colnames(ko_sum) %in% rc$SeqID]
rc.counts = rc.counts[,match(rc$SeqID, colnames(rc.counts))]

wd.counts = ko_sum[colnames(ko_sum) %in% wd$SeqID]
wd.counts = wd.counts[,match(wd$SeqID, colnames(wd.counts))]

rc.m = as.matrix(rc.counts)
wd.m = as.matrix(wd.counts)

#make deseq objects for 
dds.rc<-DESeqDataSetFromMatrix(countData=round(rc.m),colData=rc,design=~Rec_day_adj)
dds.wd<-DESeqDataSetFromMatrix(countData=round(wd.m),colData=wd,design=~Rec_day_adj)


#filter out KOs with counts of less than 20
dds.rc1 = dds.rc[rowSums(counts(dds.rc)) >= 20,]
dds.wd1 = dds.wd[rowSums(counts(dds.wd)) >= 20,]

#set D-3 (pre-ABX) as the reference level
dds.rc1$Rec_day_adj <- relevel(dds.rc1$Rec_day_adj, ref = "-3")
dds.wd1$Rec_day_adj <- relevel(dds.wd1$Rec_day_adj, ref = "-3")

#run DESeq2 on RC diet
dds.rc.out <- DESeq(dds.rc1)

rc.d3.d13 = as.data.frame(results(dds.rc.out, contrast = c("Rec_day_adj", "-13", "-3")))
rc.d3.d13_DEG = subset(rc.d3.d13, padj < 0.05) %>% #6 padj, hooray.
  subset(abs(log2FoldChange) > 2) %>% #6 at second filter
  subset(baseMean > 50) #0 at third filter
#write.csv(rc.d3.d13, "DEG_vs_preABX/rc.d13.csv")
#write.csv(rc.d3.d13_DEG, "DEG_vs_preABX/rc.d13.DEG.csv")

rc.d3.d2 = as.data.frame(results(dds.rc.out, contrast = c("Rec_day_adj", "2", "-3")))
rc.d3.d2_DEG = subset(rc.d3.d2, padj < 0.05) %>% #1816 at first filter
  subset(abs(log2FoldChange) > 2) %>% #1532 at second filter
  subset(baseMean > 50) #790 at final filter. 
#write.csv(rc.d3.d2, "DEG_vs_preABX/rc.d2.csv")
#write.csv(rc.d3.d2_DEG, "DEG_vs_preABX/rc.d2.DEG.csv")


rc.d3.d4 = as.data.frame(results(dds.rc.out, contrast = c("Rec_day_adj", "4", "-3")))
rc.d3.d4_DEG = subset(rc.d3.d4, padj < 0.05) %>% #2959 at first filter
  subset(abs(log2FoldChange) > 2) %>% #2296 at second filter
  subset(baseMean > 50) #1025 at final filter
#write.csv(rc.d3.d4, "DEG_vs_preABX/rc.d4.csv")
#write.csv(rc.d3.d4_DEG, "DEG_vs_preABX/rc.d4.DEG.csv")


rc.d3.d14 = as.data.frame(results(dds.rc.out, contrast = c("Rec_day_adj", "14", "-3")))
rc.d3.d14_DEG = subset(rc.d3.d14, padj < 0.05) %>% #224 at first filter
  subset(abs(log2FoldChange) > 2) %>% #133 at second filter
  subset(baseMean > 50) #6 at final filter
#write.csv(rc.d3.d14, "DEG_vs_preABX/rc.d14.csv")
#write.csv(rc.d3.d14_DEG, "DEG_vs_preABX/rc.d14.DEG.csv")


rc.d3.d28 = as.data.frame(results(dds.rc.out, contrast = c("Rec_day_adj", "28", "-3")))
rc.d3.d28_DEG = subset(rc.d3.d28, padj < 0.05) %>% #66 at first filter
  subset(abs(log2FoldChange) > 2) %>% #53 at second filter
  subset(baseMean > 50) #2 at final filter
#write.csv(rc.d3.d28, "DEG_vs_preABX/rc.d28.csv")
#write.csv(rc.d3.d28_DEG, "DEG_vs_preABX/rc.d28.DEG.csv")

### WD now

dds.wd.out <- DESeq(dds.wd1)

wd.d3.d13 = as.data.frame(results(dds.wd.out, contrast = c("Rec_day_adj", "-13", "-3")))
wd.d3.d13_DEG = subset(wd.d3.d13, padj < 0.05) %>% 
  subset(abs(log2FoldChange) > 2) %>% 
  subset(baseMean > 50) #115 at final filter
#write.csv(wd.d3.d13, "DEG_vs_preABX/wd.d13.csv")
#write.csv(wd.d3.d13_DEG, "DEG_vs_preABX/wd.d13.DEG.csv")

wd.d3.d2 = as.data.frame(results(dds.wd.out, contrast = c("Rec_day_adj", "2", "-3")))
wd.d3.d2_DEG = subset(wd.d3.d2, padj < 0.05) %>% 
  subset(abs(log2FoldChange) > 2) %>% 
  subset(baseMean > 50) #1110 at final filter
#write.csv(wd.d3.d2, "DEG_vs_preABX/wd.d2.csv")
#write.csv(wd.d3.d2_DEG, "DEG_vs_preABX/wd.d2.DEG.csv")

wd.d3.d14 = as.data.frame(results(dds.wd.out, contrast = c("Rec_day_adj", "14", "-3")))
wd.d3.d14_DEG = subset(wd.d3.d14, padj < 0.05) %>% 
  subset(abs(log2FoldChange) > 2) %>% 
  subset(baseMean > 50) #794 at final filter
#write.csv(wd.d3.d14, "DEG_vs_preABX/wd.d14.csv")
#write.csv(wd.d3.d14_DEG, "DEG_vs_preABX/wd.d14.DEG.csv")

wd.d3.d28 = as.data.frame(results(dds.wd.out, contrast = c("Rec_day_adj", "28", "-3")))
wd.d3.d28_DEG = subset(wd.d3.d28, padj < 0.05) %>% #1682 at first filter
  subset(abs(log2FoldChange) > 2) %>% #1277 at second filter
  subset(baseMean > 50) #338 at final filter
#write.csv(wd.d3.d28, "DEG_vs_preABX/wd.d28.csv")
#write.csv(wd.d3.d28_DEG, "DEG_vs_preABX/wd.d28.DEG.csv")

rc.d3.d13_DEG$KO = row.names()
rc.d3.d2_DEG$KO = row.names(rc.d3.d2_DEG)
rc.d3.d4_DEG$KO = row.names(rc.d3.d4_DEG)
rc.d3.d14_DEG$KO = row.names(rc.d3.d14_DEG)
rc.d3.d28_DEG$KO = row.names(rc.d3.d28_DEG)

rc.all = rbind(rc.d3.d13_DEG, rc.d3.d14_DEG, rc.d3.d2_DEG, rc.d3.d28_DEG, rc.d3.d4_DEG)
rc.degList = data.frame("KO" = rc.all$KO)
rc.degList = unique(rc.degList)
#write.csv(rc.degList, "DEG_vs_preABX/rc.DEGs.csv") #1219/1823 are unique


wd.d3.d13_DEG$KO = row.names(wd.d3.d13_DEG)
wd.d3.d2_DEG$KO = row.names(wd.d3.d2_DEG)
wd.d3.d14_DEG$KO = row.names(wd.d3.d14_DEG)
wd.d3.d28_DEG$KO = row.names(wd.d3.d28_DEG)

wd.all = rbind(wd.d3.d13_DEG, wd.d3.d14_DEG, wd.d3.d2_DEG, wd.d3.d28_DEG)
wd.degList = data.frame("KO" = wd.all$KO)
wd.degList = unique(wd.degList) #from 2357 to 1145 - lots of overlap
#write.csv(wd.degList, "DEG_vs_preABX/wd.DEGs.csv")





########## sequential day-by-day comparisons within treatment groups
rc.d2.d4 = as.data.frame(results(dds.rc.out, contrast = c("Rec_day_adj", "4", "2")))
rc.d2.d4_DEG = subset(rc.d2.d4, padj < 0.05) %>% #
  subset(abs(log2FoldChange) > 2) %>% #
  subset(baseMean > 50) #582 DEGs b/w D2 and D4
#write.csv(rc.d2.d4, "DEGs_byDay/rc.d2_4.csv")
#write.csv(rc.d2.d4_DEG, "DEGs_byDay/rc.d2_4.DEG.csv")

rc.d4.d14 = as.data.frame(results(dds.rc.out, contrast = c("Rec_day_adj", "14", "4")))
rc.d4.d14_DEG = subset(rc.d4.d14, padj < 0.05) %>% #
  subset(abs(log2FoldChange) > 2) %>% #
  subset(baseMean > 50) #986 DEGs
#write.csv(rc.d4.d14, "DEGs_byDay/rc.d4_14.csv")
#write.csv(rc.d4.d14_DEG, "DEGs_byDay/rc.d4_14.DEG.csv")

rc.d14.d28 = as.data.frame(results(dds.rc.out, contrast = c("Rec_day_adj", "28", "14")))
rc.d14.d28_DEG = subset(rc.d14.d28, padj < 0.05) %>% #
  subset(abs(log2FoldChange) > 2) %>% #
  subset(baseMean > 50) #0 DEGs


wd.d2.d14 = as.data.frame(results(dds.wd.out, contrast = c("Rec_day_adj", "14", "2")))
wd.d2.d14_DEG = subset(wd.d2.d14, padj < 0.05) %>% #
  subset(abs(log2FoldChange) > 2) %>% #
  subset(baseMean > 50) #719
#write.csv(wd.d2.d14, "DEGs_byDay/wd.d2_14.csv")
#write.csv(wd.d2.d14_DEG, "DEGs_byDay/wd.d2_14.DEG.csv")

wd.d14.d28 = as.data.frame(results(dds.wd.out, contrast = c("Rec_day_adj", "28", "14")))
wd.d14.d28_DEG = subset(wd.d14.d28, padj < 0.05) %>% #
  subset(abs(log2FoldChange) > 2) %>% #
  subset(baseMean > 50) #459
#write.csv(wd.d14.d28, "DEGs_byDay/wd.d14_28.csv")
#write.csv(wd.d14.d28_DEG, "DEGs_byDay/wd.d14_28.DEG.csv")




