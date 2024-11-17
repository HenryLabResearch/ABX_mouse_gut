library(tidyverse)
library(phyloseq)
library(microbiome)
library(pairwiseAdonis)

setwd("~/Documents/Chang-Bergelson Lab/HF Diet Abx/Sequencing/WD_16s_merged/outputs/")

otu = read.csv("wd_merged_otu.csv", row.names = 1)
meta = read.csv("wd_meta_merged_A2.csv", row.names=1) #this file has outline and shape columns for A2
tax = read.csv("wd_taxonomy.csv", row.names = 1, header=TRUE)

#filter out mitochondrial and chloroplast seqs
mito = row.names(tax)[tax$Family=="Mitochondria"]
cplast = row.names(tax)[tax$Class=="Chloroplast"]
otu.noMito = otu[!row.names(otu) %in% mito,]
tax.noMito = tax[!row.names(tax) %in% mito,]
otu.noCplast = otu.noMito[!row.names(otu.noMito) %in% cplast,]
tax.noCplast = tax.noMito[!row.names(tax.noMito) %in% cplast,]

otu_matrix = as.matrix(otu.noCplast)
tax_matrix = as.matrix(tax.noCplast)

#Use phyloseq functions to turn the imported data into phyloseq components
OTU = otu_table(otu_matrix, taxa_are_rows = TRUE)
META = sample_data(meta)
TAX = tax_table(tax_matrix)
TREE = read_tree("wd_merged_tree.nwk")
ps = phyloseq(OTU,META,TAX, TREE)

a2 = subset_samples(ps, Aim=="A2")

sample_data(a2)$Rec_day_adj = as.factor(sample_data(a2)$Rec_day_adj)
sample_data(a2)$Mouse = as.factor(sample_data(a2)$Mouse)

sample_data(a2)$Rec_day_adj = factor(sample_data(a2)$Rec_day_adj, levels=c("-31","-3","0","14","28"))
sample_data(a2)$Desc = factor(sample_data(a2)$Desc, levels=c("baseline","pre-abx","post-abx","rec-wk1","rec-wk2","rec-wk4"))

ps1 = prune_samples(sample_sums(a2)>=1000, a2)

#sum at genus level
ps3 = tax_glom(ps1, taxrank="Genus") 

#subset by timepoint
r2 = subset_samples(ps3, Desc=="rec-wk2") #D14
r4 = subset_samples(ps3, Desc=="rec-wk4") #D28

#### calculate BC distance
r2_dist = phyloseq::distance(r2, method="bray")

### ordinate on BC distance
r2_ord = ordinate(r2, "PCoA", distance = r2_dist)

#gather the PC coefficients to plot with ggplot
coefs = data.frame(r2_ord[["vectors"]])
coefs$sampleID = row.names(coefs)
meta$sampleID = row.names(meta)
coefs.meta = left_join(coefs, meta) #merge to metadata

p1 = ggplot(coefs.meta, aes(x=Axis.1, y=Axis.2,fill=Post.diet, color=Outline, shape=Shape)) +
  geom_point(size = 4 ,alpha=0.8,stroke=1) + theme_bw() +
  scale_shape_manual(values=c(0,2,17,19)) +
  scale_color_manual(values=c("#0096ff","#ff2f92","#fcada9","grey50","#77d6ff","red","blue")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none",
        #legend.title = element_text(size = 9, face="bold"), 
        #legend.text  = element_text(size = 8),
        #legend.key.size = unit(0.8, "lines"),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1)) +
  xlab("PC1 (24%)") +ylab("PC2 (17.1%)") +
  facet_wrap(Pre.diet~.)
p1

#ggsave(file="FILE.svg", plot=p1, width=7, height=3.5)




#### D28
r4_dist = phyloseq::distance(r4, method="bray")
r4_ord = ordinate(r4, "PCoA", distance = r4_dist)

coefs = data.frame(r4_ord[["vectors"]])
coefs$sampleID = row.names(coefs)
meta$sampleID = row.names(meta)
coefs.meta = left_join(coefs, meta)

p2 = ggplot(coefs.meta, aes(x=Axis.1, y=Axis.2,fill=Post.diet, color=Outline, shape=Shape)) +
  geom_point(size = 4 ,alpha=0.8,stroke=1) + theme_bw() +
  scale_shape_manual(values=c(0,2,17,19)) +
  scale_color_manual(values=c("#0096ff","#ff2f92","#fcada9","grey50","#77d6ff","red","blue")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none",
        #legend.title = element_text(size = 9, face="bold"), 
        #legend.text  = element_text(size = 8),
        #legend.key.size = unit(0.8, "lines"),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1)) +
  xlab("PC1 (21%)") +ylab("PC2 (17%)") +
  facet_wrap(Pre.diet~.)
p2

#ggsave(file="FILE.svg", plot=p2, width=7, height=3.5)










###########################################
### BETA DIVERSITY STATISTICS ###
###########################################

### Calculate mean BC distance from each treatment sample to all no-ABX PBS samples

#GENERAL OUTLINE:
#1. For each ABX sample, calculate BC distance to no-ABX controls (first compare all to RC-PBS-PBS, then to WD-PBS-PBS)
#2. Average the BC distance to no-ABX controls


#groups that change diets are WD-RC, RC-WD, WD-DIET, and RC-DIET

d14 = meta(r2)

#gather sample IDs for each treatment group
rc.wd = row.names(d14[d14$Treatment=="RC-WD",])
wd.rc = row.names(d14[d14$Treatment=="WD-RC",])
rc.rc = row.names(d14[d14$Treatment=="RC-RC",])
wd.wd = row.names(d14[d14$Treatment=="WD-WD",])
rc.micr = row.names(d14[d14$Treatment=="RC-MICR",])
wd.micr = row.names(d14[d14$Treatment=="WD-MICR",])
rc.diet = row.names(d14[d14$Treatment=="RC-DIET",])
wd.diet = row.names(d14[d14$Treatment=="WD-DIET",])
wd.abx.pbs = row.names(d14[d14$Treatment=="WD-ABX-PBS",])
rc.abx.pbs = row.names(d14[d14$Treatment=="RC-ABX-PBS",])
wd.pbs.pbs = row.names(d14[d14$Treatment=="WD-PBS-PBS",])
rc.pbs.pbs = row.names(d14[d14$Treatment=="RC-PBS-PBS",])

#take the D14 distance matrix
d14_dist = data.frame(as.matrix(phyloseq::distance(r2, method="bray")))

#compare RC-RC to RC-PBS-PBS
rc.rc_rc.pbs=NULL
for(i in 1:length(rc.rc)){
  tmp = NULL
  for(j in 1:length(rc.pbs.pbs)){
    tmp = c(tmp,d14_dist[rc.rc[i],rc.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  rc.rc_rc.pbs = c(rc.rc_rc.pbs,mean(tmp))
}
#store these in a vector, adjusting length for future merging of vectors
if(length(rc.rc_rc.pbs) < 7){
  n = 7 - length(rc.rc_rc.pbs)
  rc.rc_rc.pbs = c(rc.rc_rc.pbs, rep("NA",n))
}

#RC-WD vs RC-PBS-PBS
rc.wd_rc.pbs=NULL
for(i in 1:length(rc.wd)){
  tmp = NULL
  for(j in 1:length(rc.pbs.pbs)){
    tmp = c(tmp,d14_dist[rc.wd[i],rc.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  rc.wd_rc.pbs = c(rc.wd_rc.pbs,mean(tmp))
}
if(length(rc.wd_rc.pbs) < 7){
  n = 7 - length(rc.wd_rc.pbs)
  rc.wd_rc.pbs = c(rc.wd_rc.pbs, rep("NA",n))
}

#RC-DIET vs RC-PBS-PBS
rc.diet_rc.pbs=NULL
for(i in 1:length(rc.diet)){
  tmp = NULL
  for(j in 1:length(rc.pbs.pbs)){
    tmp = c(tmp,d14_dist[rc.diet[i],rc.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  rc.diet_rc.pbs = c(rc.diet_rc.pbs,mean(tmp))
}
if(length(rc.diet_rc.pbs) < 7){
  n = 7 - length(rc.diet_rc.pbs)
  rc.diet_rc.pbs = c(rc.diet_rc.pbs, rep("NA",n))
}

#RC-MICR vs RC-PBS-PBS
rc.micr_rc.pbs=NULL
for(i in 1:length(rc.micr)){
  tmp = NULL
  for(j in 1:length(rc.pbs.pbs)){
    tmp = c(tmp,d14_dist[rc.micr[i],rc.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  rc.micr_rc.pbs = c(rc.micr_rc.pbs,mean(tmp))
}
if(length(rc.micr_rc.pbs) < 7){
  n = 7 - length(rc.micr_rc.pbs)
  rc.micr_rc.pbs = c(rc.micr_rc.pbs, rep("NA",n))
}

#RC-ABX-PBS vs RC-PBS-PBS
rc.abx_rc.pbs=NULL
for(i in 1:length(rc.abx.pbs)){
  tmp = NULL
  for(j in 1:length(rc.pbs.pbs)){
    tmp = c(tmp,d14_dist[rc.abx.pbs[i],rc.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  rc.abx_rc.pbs = c(rc.abx_rc.pbs,mean(tmp))
}
if(length(rc.abx_rc.pbs) < 7){
  n = 7 - length(rc.abx_rc.pbs)
  rc.abx_rc.pbs = c(rc.abx_rc.pbs, rep("NA",n))
}

#WD-RC vs RC-PBS-PBS
wd.rc_rc.pbs=NULL
for(i in 1:length(wd.rc)){
  tmp = NULL
  for(j in 1:length(rc.pbs.pbs)){
    tmp = c(tmp,d14_dist[wd.rc[i],rc.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  wd.rc_rc.pbs = c(wd.rc_rc.pbs,mean(tmp))
}
if(length(wd.rc_rc.pbs) < 7){
  n = 7 - length(wd.rc_rc.pbs)
  wd.rc_rc.pbs = c(wd.rc_rc.pbs, rep("NA",n))
}

#WD-DIET vs RC-PBS-PBS
wd.diet_rc.pbs=NULL
for(i in 1:length(wd.diet)){
  tmp = NULL
  for(j in 1:length(rc.pbs.pbs)){
    tmp = c(tmp,d14_dist[wd.diet[i],rc.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  wd.diet_rc.pbs = c(wd.diet_rc.pbs,mean(tmp))
}
if(length(wd.diet_rc.pbs) < 7){
  n = 7 - length(wd.diet_rc.pbs)
  wd.diet_rc.pbs = c(wd.diet_rc.pbs, rep("NA",n))
}

#WD-MICR vs RC-PBS-PBS
wd.micr_rc.pbs=NULL
for(i in 1:length(wd.micr)){
  tmp = NULL
  for(j in 1:length(rc.pbs.pbs)){
    tmp = c(tmp,d14_dist[wd.micr[i],rc.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  wd.micr_rc.pbs = c(wd.micr_rc.pbs,mean(tmp))
}
if(length(wd.micr_rc.pbs) < 7){
  n = 7 - length(wd.micr_rc.pbs)
  wd.micr_rc.pbs = c(wd.micr_rc.pbs, rep("NA",n))
}

#WD-WD vs RC-PBS-PBS
wd.wd_rc.pbs=NULL
for(i in 1:length(wd.wd)){
  tmp = NULL
  for(j in 1:length(rc.pbs.pbs)){
    tmp = c(tmp,d14_dist[wd.wd[i],rc.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  wd.wd_rc.pbs = c(wd.wd_rc.pbs,mean(tmp))
}
if(length(wd.wd_rc.pbs) < 7){
  n = 7 - length(wd.wd_rc.pbs)
  wd.wd_rc.pbs = c(wd.wd_rc.pbs, rep("NA",n))
}

#WD-ABX-PBS vs RC-PBS-PBS
wd.abx.pbs_rc.pbs=NULL
for(i in 1:length(wd.abx.pbs)){
  tmp = NULL
  for(j in 1:length(rc.pbs.pbs)){
    tmp = c(tmp,d14_dist[wd.abx.pbs[i],rc.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  wd.abx_rc.pbs = c(wd.abx.pbs_rc.pbs,mean(tmp))
}
if(length(wd.abx.pbs_rc.pbs) < 7){
  n = 7 - length(wd.abx.pbs_rc.pbs)
  wd.abx.pbs_rc.pbs = c(wd.abx.pbs_rc.pbs, rep("NA",n))
}

#put all the comparisons into one dataframe
rc_comparisons = cbind(rc.rc_rc.pbs,rc.micr_rc.pbs,rc.abx_rc.pbs,wd.diet_rc.pbs,wd.rc_rc.pbs,
                       rc.diet_rc.pbs,rc.wd_rc.pbs,wd.micr_rc.pbs,wd.wd_rc.pbs,wd.abx_rc.pbs)

rc_comparisons = data.frame(rc_comparisons)

colnames(rc_comparisons) = c("RC-RC","RC-MICR","RC-ABX-PBS","WD-DIET","WD-RC",
                             "RC-DIET","RC-WD","WD-MICR","WD-WD","WD-ABX-PBS")

rc_comp.tidy = gather(rc_comparisons, key="Treatment",value="BC.dist")

#this file just matches up shorthand treatment groups with pre- and post-ABX treatments
trt = read.csv("trt_lookup.csv")

rc_comp.tidy = left_join(rc_comp.tidy, trt)
#write.csv(rc_comp.tidy, "a2_d14_BCdist_fromRCPBS.csv",row.names=FALSE)

#this ^ is statistically analyzed in the prism file 05.A2_BC_dist.prism


rc_comp.tidy$Treatment = factor(rc_comp.tidy$Treatment, 
                                levels=c("RC-RC","RC-MICR","RC-ABX-PBS","WD-DIET","WD-RC",
                                         "WD-WD","WD-MICR","WD-ABX-PBS","RC-DIET","RC-WD"))

#plot these data
#shows that groups that ended on RC diet are uniformly closer to no-ABX RC controls 
#than groups that ended on WD diet
ggplot(rc_comp.tidy, aes(x=Treatment, y=as.numeric(BC.dist), fill=Post.Diet)) + 
  geom_boxplot() + #ggtitle("Mean D14 BC distance from RC-PBS-PBS") +
  ylab("Mean BC distance from RC-PBS-PBS") +
  scale_fill_manual(values=c("cornflowerblue","salmon")) +
  theme(legend.title = element_text(size = 7), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.4, "lines")) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle=90),
        legend.title = element_text(face='bold'),
        #legend.position = "none",
        strip.text.x = element_text(face="bold",size=12),
        strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))




#now compare to WD-PBS-PBS
rc.rc_wd.pbs=NULL
for(i in 1:length(rc.rc)){
  tmp = NULL
  for(j in 1:length(wd.pbs.pbs)){
    tmp = c(tmp,d14_dist[rc.rc[i],wd.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  rc.rc_wd.pbs = c(rc.rc_wd.pbs,mean(tmp))
}
if(length(rc.rc_wd.pbs) < 7){
  n = 7 - length(rc.rc_wd.pbs)
  rc.rc_wd.pbs = c(rc.rc_wd.pbs, rep("NA",n))
}

rc.wd_wd.pbs=NULL
for(i in 1:length(rc.wd)){
  tmp = NULL
  for(j in 1:length(wd.pbs.pbs)){
    tmp = c(tmp,d14_dist[rc.wd[i],wd.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  rc.wd_wd.pbs = c(rc.wd_wd.pbs,mean(tmp))
}
if(length(rc.wd_wd.pbs) < 7){
  n = 7 - length(rc.wd_wd.pbs)
  rc.wd_wd.pbs = c(rc.wd_wd.pbs, rep("NA",n))
}

rc.diet_wd.pbs=NULL
for(i in 1:length(rc.diet)){
  tmp = NULL
  for(j in 1:length(wd.pbs.pbs)){
    tmp = c(tmp,d14_dist[rc.diet[i],wd.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  rc.diet_wd.pbs = c(rc.diet_wd.pbs,mean(tmp))
}
if(length(rc.diet_wd.pbs) < 7){
  n = 7 - length(rc.diet_wd.pbs)
  rc.diet_wd.pbs = c(rc.diet_wd.pbs, rep("NA",n))
}

rc.micr_wd.pbs=NULL
for(i in 1:length(rc.micr)){
  tmp = NULL
  for(j in 1:length(wd.pbs.pbs)){
    tmp = c(tmp,d14_dist[rc.micr[i],wd.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  rc.micr_wd.pbs = c(rc.micr_wd.pbs,mean(tmp))
}
if(length(rc.micr_wd.pbs) < 7){
  n = 7 - length(rc.micr_wd.pbs)
  rc.micr_wd.pbs = c(rc.micr_wd.pbs, rep("NA",n))
}

rc.abx_wd.pbs=NULL
for(i in 1:length(rc.abx.pbs)){
  tmp = NULL
  for(j in 1:length(wd.pbs.pbs)){
    tmp = c(tmp,d14_dist[rc.abx.pbs[i],wd.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  rc.abx_wd.pbs = c(rc.abx_wd.pbs,mean(tmp))
}
if(length(rc.abx_wd.pbs) < 7){
  n = 7 - length(rc.abx_wd.pbs)
  rc.abx_wd.pbs = c(rc.abx_wd.pbs, rep("NA",n))
}

wd.rc_wd.pbs=NULL
for(i in 1:length(wd.rc)){
  tmp = NULL
  for(j in 1:length(wd.pbs.pbs)){
    tmp = c(tmp,d14_dist[wd.rc[i],wd.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  wd.rc_wd.pbs = c(wd.rc_wd.pbs,mean(tmp))
}
if(length(wd.rc_wd.pbs) < 7){
  n = 7 - length(wd.rc_wd.pbs)
  wd.rc_wd.pbs = c(wd.rc_wd.pbs, rep("NA",n))
}

wd.diet_wd.pbs=NULL
for(i in 1:length(wd.diet)){
  tmp = NULL
  for(j in 1:length(wd.pbs.pbs)){
    tmp = c(tmp,d14_dist[wd.diet[i],wd.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  wd.diet_wd.pbs = c(wd.diet_wd.pbs,mean(tmp))
}
if(length(wd.diet_wd.pbs) < 7){
  n = 7 - length(wd.diet_wd.pbs)
  wd.diet_wd.pbs = c(wd.diet_wd.pbs, rep("NA",n))
}

wd.micr_wd.pbs=NULL
for(i in 1:length(wd.micr)){
  tmp = NULL
  for(j in 1:length(wd.pbs.pbs)){
    tmp = c(tmp,d14_dist[wd.micr[i],wd.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  wd.micr_wd.pbs = c(wd.micr_wd.pbs,mean(tmp))
}
if(length(wd.micr_wd.pbs) < 7){
  n = 7 - length(wd.micr_wd.pbs)
  wd.micr_wd.pbs = c(wd.micr_wd.pbs, rep("NA",n))
}

wd.wd_wd.pbs=NULL
for(i in 1:length(wd.wd)){
  tmp = NULL
  for(j in 1:length(wd.pbs.pbs)){
    tmp = c(tmp,d14_dist[wd.wd[i],wd.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  wd.wd_wd.pbs = c(wd.wd_wd.pbs,mean(tmp))
}
if(length(wd.wd_wd.pbs) < 7){
  n = 7 - length(wd.wd_wd.pbs)
  wd.wd_wd.pbs = c(wd.wd_wd.pbs, rep("NA",n))
}

wd.abx.pbs_wd.pbs=NULL
for(i in 1:length(wd.abx.pbs)){
  tmp = NULL
  for(j in 1:length(wd.pbs.pbs)){
    tmp = c(tmp,d14_dist[wd.abx.pbs[i],wd.pbs.pbs[j]]) #this collects the distance between this treatment group and all the PBS control dots
  }
  wd.abx_wd.pbs = c(wd.abx.pbs_wd.pbs,mean(tmp))
}
if(length(wd.abx.pbs_wd.pbs) < 7){
  n = 7 - length(wd.abx.pbs_wd.pbs)
  wd.abx.pbs_wd.pbs = c(wd.abx.pbs_wd.pbs, rep("NA",n))
}

wd_comparisons = cbind(rc.rc_wd.pbs,rc.micr_wd.pbs,rc.abx_wd.pbs,wd.diet_wd.pbs,wd.rc_wd.pbs,
                       rc.diet_wd.pbs,rc.wd_wd.pbs,wd.micr_wd.pbs,wd.wd_wd.pbs,wd.abx_wd.pbs)

wd_comparisons = data.frame(wd_comparisons)

colnames(wd_comparisons) = c("RC-RC","RC-MICR","RC-ABX-PBS","WD-DIET","WD-RC",
                             "RC-DIET","RC-WD","WD-MICR","WD-WD","WD-ABX-PBS")

wd_comp.tidy = gather(wd_comparisons, key="Treatment",value="BC.dist")

trt = read.csv("trt_lookup.csv")

wd_comp.tidy = left_join(wd_comp.tidy, trt)
#write.csv(wd_comp.tidy, "a2_d14_BCdist_fromWDPBS.csv",row.names=FALSE)

# ^this is statistically analyzed in 05.A2_BC_dist.prism

wd_comp.tidy$Treatment = factor(wd_comp.tidy$Treatment, 
                                   levels=c("RC-RC","RC-MICR","RC-ABX-PBS","WD-DIET","WD-RC",
                                            "WD-WD","WD-MICR","WD-ABX-PBS","RC-DIET","RC-WD"))


p2 = ggplot(wd_comp.tidy, aes(x=Treatment, y=as.numeric(BC.dist), fill=Post.Diet)) + 
  geom_boxplot() + #ggtitle("Mean D14 BC distance from WD-PBS-PBS") +
  ylab("Mean BC distance from WD-PBS-PBS") + ylim(c(0.6,1)) +
  scale_fill_manual(values=c("cornflowerblue","salmon")) +
  theme(legend.title = element_text(size = 7), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.4, "lines")) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle=90),
        legend.title = element_text(face='bold'),
        #legend.position = "none",
        strip.text.x = element_text(face="bold",size=12),
        strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p2




#################################################
### PERMANOVA: 
### BC dist ~ Post.diet*ABX*Post.micr + Pre-diet
#################################################

#we're including an interaction effect for post-diet and antibiotics bc it seems like 
#the effect of antibiotics is different depending on post-diet
r2_meta = meta(r2)
anova(betadisper(r2_dist, r2_meta$Post.micr)) #not diff for pre.diet, post.diet, or post.micr
out = adonis2(r2_dist ~ r2_meta$Post.diet*r2_meta$Abx*r2_meta$Post.micr + r2_meta$Pre.diet, data = r2_meta, perm=9999)
out
posthoc = pairwise.adonis(r2_dist, r2_meta$Treatment, p.adjust.m="BH")
posthoc


