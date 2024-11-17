library(tidyverse)
library(phyloseq)
library(microbiome)
library(picante)

setwd("/DIRECTORY/")

otu = read.csv("wd_merged_otu.csv", row.names = 1)
meta = read.csv("wd_meta_merged.csv", row.names=1)
tax = read.csv("wd_taxonomy.csv", row.names = 1, header=TRUE)

otu_matrix = as.matrix(otu)
tax_matrix = as.matrix(tax)

# Use phyloseq functions to turn the imported data into phyloseq components
OTU = otu_table(otu_matrix, taxa_are_rows = TRUE)
META = sample_data(meta)
TAX = tax_table(tax_matrix)
TREE = read_tree("wd_merged_tree.nwk")
ps = phyloseq(OTU,META,TAX, TREE)

sample_data(ps)$Rec_day_adj = as.factor(sample_data(ps)$Rec_day_adj)
sample_data(ps)$Mouse = as.factor(sample_data(ps)$Mouse)

ps1 = prune_samples(sample_sums(ps)>=1000, ps)
a2 = subset_samples(ps1, Aim=="A2")

shannon <-microbiome::alpha(a2, index = "shannon")
richness = microbiome::richness(a2)
pd = picante::pd(t(otu_table(a2)), TREE)
pd$sampleID = row.names(pd)

dm = cbind(shannon, richness, pd)
dm$sampleID= rownames(dm)

meta$sampleID = row.names(meta)
dm1 = left_join(dm, meta, by="sampleID")

#to check how many Ns per group/timepoint
dm1 %>% group_by(Treatment, Rec_day_adj) %>% tally() %>% print(n=66)

#write.csv(dm1, "a2_dm.csv",row.names=FALSE)


#Figure 4C, Figure S7C, and statistical analyses were generated in prism


















#let's just look through D14
dm2 = subset(dm1, Rec_day < 16)


a2.summ.h = dm2 %>% group_by(Rec_day_adj, Treatment, Cohort) %>% 
  summarise(
    mean_H = mean(diversity_shannon),
    sd = sd(diversity_shannon, na.rm=TRUE))


ggplot(a2.summ.h) + aes(x = Rec_day_adj, y = mean_H, group = Treatment, col = Treatment) +
  geom_line() + geom_point() + 
  ylab("Shannon Diversity") + 
  facet_wrap(Cohort~.)


######### richness now
a2.summ.rich = dm2 %>% group_by(Rec_day_adj, Treatment, Cohort, Pre.diet, Post.diet, Post.micr, Abx) %>% 
  summarise(
    mean_taxa = mean(observed),
    sd = sd(observed, na.rm=TRUE))

ggplot(a2.summ.rich) + aes(x = Rec_day_adj, y = mean_taxa, group = Treatment, col = Treatment) +
  geom_line() + geom_point() + 
  ylab("Observed ASVs") + 
  facet_wrap(Cohort~.)


ggplot(a2.summ.rich) + aes(x = Rec_day_adj, y = mean_taxa, group = Treatment, col = Treatment) +
  geom_smooth() + geom_point() + 
  ylab("Observed ASVs") + 
  facet_wrap(Treatment~.)

ggplot(a2.summ.rich) + aes(x = Rec_day_adj, y = mean_taxa, group = Treatment, col = Treatment) +
  geom_smooth() + geom_point() + 
  ylab("Observed ASVs") + 
  facet_wrap(Post.diet~.)

#write.csv(a2.summ.rich, "a2_richness.csv",row.names=FALSE)
















a2.summ.rich$trt.coh = paste0(a2.summ.rich$Treatment, a2.summ.rich$Cohort)

###plot all on the same graph
p2 = ggplot(a2.summ.rich) + aes(x = Rec_day_adj, y = mean_taxa, group = trt.coh, col = Treatment, linetype=Cohort) +
  geom_line() + geom_point() + 
  #geom_errorbar(aes(ymin=mean_H-sd, ymax=mean_H+sd),
  #              #position=position_dodge(0.5),
  #              alpha=0.6, width=3) +
  ylab("Observed ASVs") + 
  scale_x_continuous(breaks=c(-31,-13,-3,0,3,7,14,21, 28, 35, 49, 63)) +
  scale_y_continuous(breaks=c(50,100,150,200,250,300)) +
  scale_color_manual(values=c("blue","turquoise","red3","salmon")) +
  scale_linetype_manual(values=c(1,2,3)) +
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
#ggsave(file="A1_richness_byCohort.svg", plot=p2, width=6, height=4)


a1.sub.summ = a1.sub %>% group_by(Rec_day_adj, Treatment) %>% 
  summarise(
    mean_H = mean(diversity_shannon),
    sd = sd(diversity_shannon, na.rm=TRUE))

ggplot(a1.sub.summ) + aes(x = Rec_day_adj, y = mean_H, group = Treatment, col = Treatment) +
  geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=mean_H-sd, ymax=mean_H+sd),
                position=position_dodge(0.05),
                alpha=0.4) +
  ylab("Shannon Diversity") + 
  scale_x_continuous(breaks=c(-13,-3, 0,3,7, 11, 14, 16, 18, 21, 28))








ggplot(a1.summ) + aes(x = Rec_day_adj, y = mean_H, group = Treatment, col = Treatment) +
  geom_line(linewidth=.8) + 
  geom_point() +
  geom_ribbon(aes(ymin=mean_H-sd, ymax=mean_H+sd, fill=Treatment), alpha = 0.25, col=NA, show.legend = FALSE) +
  ylab("Shannon Diversity") + xlab("Day") +
  facet_grid(Cohort~.) +
  scale_x_continuous(limits = c(-31, 28), breaks=c(-31,-3,0,3,5,6,7,8,11,14, 18, 21, 28, 36, 50, 64))

dm2.summ = dm1 %>% group_by(Rec_day_adj, Treatment, Aim) %>% 
  summarise(
    mean_H = mean(diversity_shannon),
    sd = sd(diversity_shannon, na.rm=TRUE))

a1.summ2 = subset(dm2.summ, Aim=="A1")
a2.summ2 = subset(dm2.summ, Aim=="A2")

#I think i need to adjust the metadata so the days match
ggplot(a1.summ2) + aes(x = Rec_day_adj, y = mean_H, group = Treatment, col = Treatment) +
  geom_line(size=.8) + 
  geom_point() +
  geom_ribbon(aes(ymin=mean_H-sd, ymax=mean_H+sd, fill=Treatment), alpha = 0.25, col=NA, show.legend = FALSE) +
  ylab("Shannon Diversity") + xlab("Day") +
  scale_x_continuous(limits = c(-31, 28), breaks=c(-31,-3,0,3,5,8,11,15, 18, 21, 28, 36, 50, 64))



dm1_richness.summ = dm1 %>% group_by(Rec_day_adj, Treatment, Aim, Cohort) %>% 
  summarise(
    mean_taxa = mean(observed),
    sd = sd(observed, na.rm=TRUE))

ggplot(dm1_richness.summ) + aes(x = Rec_day_adj, y = mean_taxa, group = Treatment, col = Treatment) +
  geom_ribbon(aes(ymin=mean_taxa-sd, ymax=mean_taxa+sd, fill=Treatment), alpha = 0.15, col=NA, show.legend = FALSE) +
  geom_line(size=0.8) + 
  #scale_color_manual(values=cols) + scale_fill_manual(values=cols) +
  geom_point() + facet_grid(Cohort~Aim) +
  ylab("Observed Taxa") + xlab("Day") +
  scale_x_continuous(breaks=c(-13,-3,0,3,5,8,11,15, 18, 22, 29, 36, 50, 64))
  

ggplot(dm1) + aes(x = day, y = chao1, col = treatment) +
  geom_point() + geom_smooth() +
  ylab("Chao1") + xlab("Day") +
  labs(col="Diet") + scale_x_continuous(breaks=c(-13,-3,0,3,5,8,11,15)) +
  ggtitle("Richness")

ggplot(dm1) + aes(x = new_day, y = observed, col = treatment) +
  geom_point() + geom_smooth() +
  ylab("Observed Taxa") + xlab("Day") +
  labs(col="Diet") + scale_x_continuous(breaks=c(-13,-3,0,3,5,8,11,15)) +
  ggtitle("Richness")



ggplot(a1) + aes(x = Rec_day_adj, y = PD, col = Treatment) +
  geom_point() + geom_smooth() +
  ylab("Shannon Diversity") + xlab("Day") +
  labs(col="Diet") + scale_x_continuous(breaks=c(-13,-3,0,3,5,8,11,15, 18, 22, 29, 36, 50, 64)) +
  ggtitle("Faith's PD") + ylim(c(0,75)) +
  facet_grid(Cohort~.)

a1_pd.summ = a1 %>% group_by(Rec_day_adj, Treatment, Cohort) %>% 
  summarise(
    mean_pd = mean(PD),
    sd = sd(PD, na.rm=TRUE))

ggplot(a1_pd.summ) + aes(x = Rec_day_adj, y = mean_pd, group = Treatment, col = Treatment) +
  geom_ribbon(aes(ymin=mean_pd-sd, ymax=mean_pd+sd, fill=Treatment), alpha = 0.15, col=NA, show.legend = FALSE) +
  geom_line(size=0.8) + 
  #scale_color_manual(values=cols) + scale_fill_manual(values=cols) +
  geom_point() + facet_grid(Cohort~.) +
  ylab("Faith's PD") + xlab("Day") +
  scale_x_continuous(breaks=c(-13,-3,0,3,5,8,11,15, 18, 22, 29, 36, 50, 64))


#I think for the final figure, we're just going to show richness (others can go in supp)
a1.c13 = subset(a1, Cohort!="C4") %>% subset(Cohort !="C5") %>% subset(Cohort != "C6")

a1c13_richness.summ = a1.c13 %>% group_by(Rec_day_adj, Treatment, Cohort) %>% 
  summarise(
    mean_taxa = mean(observed),
    sd = sd(observed, na.rm=TRUE))

#I might need to make a treatment-cohort variable
a1c13_richness.summ$trt.coh = paste0(a1c13_richness.summ$Treatment, a1c13_richness.summ$Cohort)


ggplot(a1c13_richness.summ) + aes(x = Rec_day_adj, y = mean_taxa, group = trt.coh, col = Treatment, lty=Cohort) +
  #geom_ribbon(aes(ymin=mean_taxa-sd, ymax=mean_taxa+sd, fill=Treatment), alpha = 0.15, col=NA, show.legend = FALSE) +
  geom_line(linewidth=0.8) + 
  #scale_color_manual(values=cols) + scale_fill_manual(values=cols) +
  geom_point() + 
  ylab("Observed Taxa") + xlab("Day") +
  scale_x_continuous(breaks=c(-13,-3,0,3,5,8,11,15, 18, 22, 29, 36, 50, 64))

#figure out which sample is the cohort 2 RC-ABX D25 one that's causing that dip
#and why there's the peak in C2 RC-PBS D14

#tbh, consider doing richness with rarefied data?

