library(tidyverse)
library(phyloseq)
library(microbiome)
library(picante)

setwd("/DIRECTORY/")

#STEP 1: calculate diversity metrics for both aims

#import data
otu = read.csv("wd_merged_otu.csv", row.names = 1)
meta = read.csv("wd_meta_merged.csv", row.names=1)
tax = read.csv("wd_taxonomy.csv", row.names = 1, header=TRUE)

otu_matrix = as.matrix(otu)
tax_matrix = as.matrix(tax)

#convert to phyloseq objects and import to phyloseq
OTU = otu_table(otu_matrix, taxa_are_rows = TRUE)
META = sample_data(meta)
TAX = tax_table(tax_matrix)
TREE = read_tree("wd_merged_tree.nwk")
ps = phyloseq(OTU,META,TAX,TREE)

#Rec_day_adj = recovery day = days post-abx (detailed description in ReadMe.txt)
#adjust numeric values to factors for Rec Day and mouse ID
sample_data(ps)$Rec_day_adj = as.factor(sample_data(ps)$Rec_day_adj)
sample_data(ps)$Mouse = as.factor(sample_data(ps)$Mouse)

ps1 = subset_samples(ps, nseqs > 5000)

#subset the samples for the initial WD Resilience experiments (aim 1)
a1 = subset_samples(ps1, Aim=="A1")

#Cohorts 4-6 were excluded due to experimental errors :( details in ReadMe.txt
a1c123 = subset_samples(a1, Cohort!= "C4")
a1c123 = subset_samples(a1c123, Cohort!= "C5")
a1c123 = subset_samples(a1c123, Cohort!= "C6")

#calculate diversity metrics
shannon.123 <-microbiome::alpha(a1c123, index = "shannon")
richness.123 = microbiome::richness(a1c123)
evenness.123 = microbiome::evenness(a1c123, "all")
pd.123 = picante::pd(t(otu_table(a1c123)), TREE)
pd.123$sampleID = row.names(pd.123)

#combine the diversity metrics into one dataframe
dm.123 = cbind(shannon.123, richness.123, evenness.123)
dm.123$sampleID= rownames(dm.123)
dm.123 = left_join(dm.123, pd.123)

meta$sampleID = row.names(meta)

#merge the diversity metrics with the metadata for future reference
dm.123 = left_join(dm.123, meta, by="sampleID")

###subset only the days that are common to all cohorts
good.days = c(-3,0,1,2,3,4,5,7,8,11,14,18,21,28)

a1.sub = subset(dm.123, Rec_day_adj %in% good.days)

a1.sub %>% group_by(Treatment, Rec_day_adj) %>% tally() %>% print(n=40)
#n=4-14

#write.csv(a1.sub, "a1_dm.csv")

#the main text figures were generated in Prism using this^ output


#Figures S1D-F were generated below:
#FIRST PLOTS: cohort-split Shannon, richness, and PD

#calculate mean Shannon for each timepoint, treatment group, and cohort
a1.summ.h = a1.sub %>% group_by(Rec_day_adj, Treatment, Cohort) %>% 
  summarise(
    mean_H = mean(diversity_shannon),
    sd = sd(diversity_shannon, na.rm=TRUE))

a1.summ.h$trt.coh = paste0(a1.summ.h$Treatment, a1.summ.h$Cohort)

###plot all cohorts faceted on the same graph
p1=ggplot(a1.summ.h) + aes(x = Rec_day_adj, y = mean_H, group = trt.coh, col = Treatment) +
  geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=mean_H-sd, ymax=mean_H+sd),
                #position=position_dodge(0.5),
                alpha=0.6, width=3) +
  ylab("Shannon Diversity") + 
  facet_grid(Cohort~.) +
  scale_x_continuous(breaks=c(-31,-13,-3,0,3,7,14,21, 28, 35, 49, 63)) +
  scale_color_manual(values=c("blue","turquoise","red3","salmon")) +
  #scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        axis.text.y = element_text(size=12, face="bold", color='black'),
        #axis.text.x = element_blank(),
        legend.title = element_text(face='bold'))
p1
#ggsave(file="FILE.svg", plot=p1, width=7, height=5)



######### richness now
a1.summ.rich = a1.sub %>% group_by(Rec_day_adj, Treatment, Cohort) %>% 
  summarise(
    mean_taxa = mean(observed),
    sd = sd(observed, na.rm=TRUE))

a1.summ.rich$trt.coh = paste0(a1.summ.rich$Treatment, a1.summ.rich$Cohort)

###plot all on the same graph
p2 = ggplot(a1.summ.rich) + aes(x = Rec_day_adj, y = mean_taxa, group = trt.coh, col = Treatment) +
  geom_line() + geom_point() + 
  facet_grid(Cohort~.) +
  geom_errorbar(aes(ymin=mean_taxa-sd, ymax=mean_taxa+sd),
                #position=position_dodge(1),
                alpha=0.6, width=3) +
  ylab("Observed ASVs") + 
  scale_x_continuous(breaks=c(-31,-13,-3,0,3,7,14,21, 28, 35, 49, 63)) +
  scale_y_continuous(breaks=c(50,100,150,200,250,300)) +
  scale_color_manual(values=c("blue","turquoise","red3","salmon")) +
  #scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        axis.text.y = element_text(size=12, face="bold", color='black'),
        #axis.text.x = element_blank(),
        legend.title = element_text(face='bold'))
p2
#ggsave(file="FILE.svg", plot=p2, width=7, height=5)


######### PD now
a1.summ.pd = a1.sub %>% group_by(Rec_day_adj, Treatment, Cohort) %>% 
  summarise(
    mean_taxa = mean(PD),
    sd = sd(PD, na.rm=TRUE))

a1.summ.pd$trt.coh = paste0(a1.summ.pd$Treatment, a1.summ.pd$Cohort)

###plot all on the same graph
p3 = ggplot(a1.summ.pd) + aes(x = Rec_day_adj, y = mean_taxa, group = trt.coh, col = Treatment) +
  geom_line() + geom_point() + 
  facet_grid(Cohort~.) +
  geom_errorbar(aes(ymin=mean_taxa-sd, ymax=mean_taxa+sd),
                #position=position_dodge(1),
                alpha=0.6, width=3) +
  ylab("Faith's PD") + 
  scale_x_continuous(breaks=c(-31,-13,-3,0,3,7,14,21, 28, 35, 49, 63)) +
  #scale_y_continuous(breaks=c(50,100,150,200,250,300)) +
  scale_color_manual(values=c("blue","turquoise","red3","salmon")) +
  #scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        axis.text.y = element_text(size=12, face="bold", color='black'),
        #axis.text.x = element_blank(),
        legend.title = element_text(face='bold'))
p3
#ggsave(file="FILE.svg", plot=p3, width=7, height=5)





