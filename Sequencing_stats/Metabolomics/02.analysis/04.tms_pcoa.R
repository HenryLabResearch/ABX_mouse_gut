library(tidyverse)
library(vegan)
library(microbiome)
library(pairwiseAdonis)
library(phyloseq)

setwd("DIRECTORY/")

#TMS panel
gcms = read.csv("TMS.csv", check.names = FALSE, row.names=1)
#set the NAs to 0
gcms[is.na(gcms)] = 0

meta = read.csv("metabolomics_meta.csv",row.names=1)

# Use phyloseq functions to turn the imported data into phyloseq components
OTU = otu_table(gcms, taxa_are_rows = FALSE)
META = sample_data(meta)
ps = phyloseq(OTU,META)

#factor the numeric values
sample_data(ps)$Rec_day_adj = factor(sample_data(ps)$Rec_day_adj, levels=c("-3","0","3","5","7","11","14","28"))
sample_data(ps)$Cage = as.factor(sample_data(ps)$Cage)
sample_data(ps)$Mouse = as.factor(sample_data(ps)$Mouse)

#first analyze fecal samples
fecal = subset_samples(ps, Type=="fecal")

#calculate pairwise BC distance
fecal_dist = phyloseq::distance(fecal, method="bray")
#ordinate on BC distance
fecal_ord = ordinate(fecal, "PCoA", distance = fecal_dist)

#Figure S3B    
plot_ordination(fecal, fecal_ord, color = "Rec_day_adj",shape="Treatment") +
    geom_point(size = 4, alpha=0.7) + theme_bw() +
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(size = 1),
      panel.background = element_blank()) +
  facet_wrap(Treatment~.) + ggtitle("TMS metabolomics panel; fecal") +
  theme(plot.title = element_text(face="bold"),
          axis.title.x = element_text(face="bold"),
          axis.title.y = element_text(face="bold"),
          axis.text.y = element_text(size=12, face="bold", color='black'),
          axis.text.x = element_text(size = 12, face='bold',color='black'),
          legend.title = element_text(face='bold'))


#PERMANOVA clustering analysis
#Table S5B
fecal_meta = meta(fecal)
anova(betadisper(fecal_dist, fecal_meta$Treatment)) #sig. diff
out = adonis2(fecal_dist ~ fecal_meta$Rec_day_adj*fecal_meta$Treatment, data = fecal_meta, perm=9999)
out

#posthoc
rc = subset_samples(fecal, Treatment=="RC-ABX")
wd = subset_samples(fecal, Treatment=="WD-ABX")

rc_dist = phyloseq::distance(rc, method="bray")
wd_dist = phyloseq::distance(wd, method="bray")

rc_meta = meta(rc)
wd_meta = meta(wd)

rc_meta$day_diet = paste0(rc_meta$Treatment, "_", rc_meta$Rec_day_adj)
rc.posthoc = pairwise.adonis(rc_dist, rc_meta$day_diet, p.adjust.m="BH")
rc.posthoc
#write.csv(rc.posthoc, "TMS_posthoc_RC.csv")

wd_meta$day_diet = paste0(wd_meta$Treatment, "_", wd_meta$Rec_day_adj)
wd.posthoc = pairwise.adonis(wd_dist, wd_meta$day_diet, p.adjust.m="BH")
wd.posthoc
#write.csv(wd.posthoc, "TMS_posthoc_WD.csv")



### now cecal samples 

cecal = subset_samples(ps, Type=="cecal")
cecal_dist = phyloseq::distance(cecal, method="bray")
    
cecal_ord = ordinate(cecal, "PCoA", distance = cecal_dist)


#plotted by extracting coefficiencts and using ggplot for better flexibility
#with shapes and formatting
coefs = data.frame(cecal_ord[["vectors"]])
coefs$sampleID = row.names(coefs)
meta$sampleID = row.names(meta)
coefs.meta = left_join(coefs, meta)

coefs.meta$Rec_day_adj = factor(coefs.meta$Rec_day_adj, levels=c("-3","14","28"))

#Figure S3C
ggplot(coefs.meta, aes(x=Axis.1, y=Axis.2,color=Rec_day_adj,shape=Treatment)) +
  geom_point(size = 4) + theme_bw() +
  scale_shape_manual(values=c(16,1,17,2)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        panel.background = element_blank()) +
  facet_wrap(Treatment~.) + ggtitle("TMS metabolomics panel; cecal") +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        axis.text.y = element_text(size=12, face="bold", color='black'),
        axis.text.x = element_text(size = 12, face='bold',color='black'),
        legend.title = element_text(face='bold'))
  