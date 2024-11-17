#goal of this script is to take specific subsets of genes and to plot their 
#relative abundances over time

library(tidyverse)
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(gcookbook)
library(scales)
library(RColorBrewer)

#you have to do this at the KO level, bc there are too many other non-DEG KOs within each kcat

setwd("DIRECTORY/")

#import counts per KO
counts = read.csv("counts_ko_sum.csv")

#convert to relative abundances
rel.cov = counts
for(i in 2:ncol(rel.cov)){
  tot = sum(rel.cov[,i])
  rel.cov[,i] = rel.cov[,i]/tot
}

row.names(rel.cov) = rel.cov$KO
rel.cov1 = rel.cov[,-1] #get rid of the column for KO

ko.m = as.matrix(rel.cov1)

#read in sample metadata
meta = read.csv("meta.csv")
row.names(meta) = meta$SeqID

OTU = otu_table(ko.m, taxa_are_rows = TRUE)
META = sample_data(meta)
ps = phyloseq(OTU,META)

sample_data(ps)$Rec_day_adj = factor(sample_data(ps)$Rec_day_adj, levels=c("-13","-3","2","4","14","28"))

#subset samples present at all timepoints
complete = c(1507, 1512, 1946, 1953, 5733, 5735)
ps.cpt = subset_samples(ps, Mouse %in% complete)

#extract counts/metadata in tidy form
ps_df = psmelt(ps.cpt)
colnames(ps_df)[colnames(ps_df)=="OTU"] = "KO"

#filter out the pre-diet-acclimation timepoint
ps_df1 = subset(ps_df, Rec_day_adj != "-13")

#this file converts between KOs and ECs
ec = read.csv("ko.ec.csv")



###############################################################
### 1. alpha-galactosidase genes
###############################################################


#this file includes all the ECs that map to agal metabolism in the dbCAN database
agal = read.csv("agal_genes.csv")

#match up KOs in df with ECs
ps2 = left_join(ps_df1, ec, multiple="all")

#filter to only include KOs that map to agal metabolism
agal.hits = ps2[ps2$EC %in% agal$EC,]

#average KO abundances across mice
agal.mean = agal.hits %>% group_by(Rec_day_adj,Treatment, KO) %>%
  summarize(meanAbun = mean(Abundance))

#sum up all the agal genes to simplify
agal.sum = agal.hits %>% group_by(Rec_day_adj, Treatment, Mouse) %>%
  summarize(sumAbun = sum(Abundance))

ggplot(agal.sum, aes(x=Rec_day_adj, y=sumAbun, group=Treatment,color=Treatment)) +
  geom_point() + geom_smooth() +
  #facet_wrap(Treatment~.,nrow=2) +
  scale_y_log10()
#shell yeah man

#then go plot this in prism
#write.csv(agal.sum, "agal.sum.csv",row.names=FALSE)



###################################################################
### 2. polysaccharide genes
### manually picked 4 polysaccharide substrates that were in our 
### metabolomics dataset: starch, cellulose, arabinan, and xylose
###################################################################

psac = read.csv("psacchs.csv")

#pick out KOs that map to metabolism of these polysaccharides
psac.hits = ps2[ps2$EC %in% psac$EC,]

#this just merges to include the column that says which susbtrate is metabolized
#by which genes in this dataset
psac.hits1 = left_join(psac.hits, psac, multiple="all")

#average across mice and sum up the genes for each substrate
psac.sum = psac.hits1 %>% group_by(Rec_day_adj, Treatment, Mouse, Substrate) %>%
  summarize(sumAbun = sum(Abundance))

ggplot(psac.sum, aes(x=Rec_day_adj, y=sumAbun, group=Substrate,color=Substrate)) +
  geom_point() + geom_smooth() +
  facet_wrap(Treatment~.,nrow=2) +
  scale_y_log10()

#write.csv(psac.sum, "psac.hits.csv",row.names=FALSE)
#plot this in prism



#######################################################################
### 3. SCFA genes
### specifically, looking at acetate, butyrate, propionate, succinate
#######################################################################

scfas = read.csv("SCFA_genes.csv")

unique(scfas$COMPOUND)

ace = subset(scfas, COMPOUND=="Acetate")
ace.hits = ps_df1[ps_df1$KO %in% ace$KO,]

ace.sum = ace.hits %>% group_by(Rec_day_adj, Treatment, Mouse) %>%
  summarize(sumAbun = sum(Abundance))

ggplot(ace.sum, aes(x=Rec_day_adj, y=sumAbun, group=Treatment,color=Treatment)) +
  geom_point() + geom_smooth() +
  #facet_wrap(Treatment~.,nrow=2) +
  scale_y_log10()


but = subset(scfas, COMPOUND=="Butyrate")
but.hits = ps_df1[ps_df1$KO %in% but$KO,]

but.sum = but.hits %>% group_by(Rec_day_adj, Treatment, Mouse) %>%
  summarize(sumAbun = sum(Abundance))

ggplot(but.sum, aes(x=Rec_day_adj, y=sumAbun, group=Treatment,color=Treatment)) +
  geom_point() + geom_smooth() #+
  #facet_wrap(Treatment~.,nrow=2) +
  #scale_y_log10()

prop = subset(scfas, COMPOUND=="propionate")
prop.hits = ps_df1[ps_df1$KO %in% prop$KO,]

prop.sum = prop.hits %>% group_by(Rec_day_adj, Treatment, Mouse) %>%
  summarize(sumAbun = sum(Abundance))

ggplot(prop.sum, aes(x=Rec_day_adj, y=sumAbun, group=Treatment,color=Treatment)) +
  geom_point() + geom_smooth() #+
#facet_wrap(Treatment~.,nrow=2) +
#scale_y_log10()

#merge and write the scfas
ace.sum$compound = "acetate"
but.sum$compound = "butyrate"
prop.sum$compound = "propionate"

scfa.hits.all = rbind(ace.sum, but.sum, prop.sum)
#write.csv(scfa.hits.all, "scfa.hits.csv",row.names=FALSE)
#analyze this in prism




#######################################################################
### 3. Bile acid metabolism genes:
### K00076 = hdhA; 7-alpha-hydroxysteroid dehydrogenase [EC:1.1.1.159]
### K01442 = cbh; choloylglycine hydrolase [EC:3.5.1.24] aka BSH
### K15868 = baiB; bile acid-coenzyme A ligase [EC:6.2.1.7]
#######################################################################

#these KOs map to HDHA, BaiB, and BSH
ba = c("K00076","K01442","K15868")
ba.hits = subset(ps_df1, KO %in% ba)

ba.mean = ba.hits %>% group_by(Rec_day_adj,Treatment, KO) %>%
  summarize(meanAbun = mean(Abundance))

ba.mean$ko.trt = paste(ba.mean$KO, ba.mean$Treatment, sep="_")

ggplot(ba.mean, aes(x=Rec_day_adj, y=meanAbun, group=ko.trt, color=Treatment,linetype=KO)) +
  geom_point() + geom_line() +
  #facet_wrap(Treatment~.,nrow=2) +
  scale_y_log10()

#write.csv(ba.hits, "ba.hits.csv",row.names=FALSE)

