library(tidyverse)
library(phyloseq)
library(microbiome)
library(pairwiseAdonis)

setwd("/DIRECTORY/")


####import the data
otu = read.csv("wd_merged_otu.csv", header=TRUE, row.names = 1)
meta = read.csv("wd_meta_merged.csv", header=TRUE)
row.names(meta) = meta$sampleID
tax = read.csv("wd_taxonomy.csv", row.names = 1, header=TRUE)

#filter mitochondrial seqs
mito = row.names(tax)[tax$Family=="Mitochondria"]
otu.noMito = otu[!row.names(otu) %in% mito,]
tax.noMito = tax[!row.names(tax) %in% mito,]

#convert to matrix
otu_matrix = as.matrix(otu.noMito)
tax_matrix = as.matrix(tax.noMito)

#convert to phyloseq objects
OTU = otu_table(otu_matrix, taxa_are_rows = TRUE)
META = sample_data(meta)
TAX = tax_table(tax_matrix)
TREE = read_tree("wd_merged_tree.nwk")

ps = phyloseq(OTU,META,TAX,TREE)

#subset the initial WD-resilience experiments
a1 = subset_samples(ps, Aim=="A1")

#cohorts 4-6 excluded due to experimental error
a1c123 = subset_samples(a1, Cohort!= "C4")
a1c123 = subset_samples(a1c123, Cohort!= "C5")
a1c123 = subset_samples(a1c123, Cohort!= "C6")

#get rid of samples with <5000 reads
ps1 = prune_samples(sample_sums(a1c123)>=5000, a1c123)
#get rid of taxa with low counts
ps2 = prune_taxa(taxa_sums(ps1) > 4, ps1)
#convert recovery day to factor
sample_data(ps2)$Rec_day_adj = as.factor(sample_data(ps2)$Rec_day_adj)

#aggregate at genus level
ps_gen = tax_glom(ps2, taxrank="Genus") 

#calculate BC distance
ps_gen_dist = phyloseq::distance(ps_gen, method="bray")
#ordinate on BC distance
ps_gen_ord = ordinate(ps_gen, "PCoA", distance= ps_gen_dist)

#gather the PC coefficients from the ordination 
coefs = data.frame(ps_gen_ord[["vectors"]])
coefs$sampleID = row.names(coefs)
coefs.meta = left_join(coefs, meta)

#set recovery day as factor
coefs.meta$Rec_day_adj = factor(coefs.meta$Rec_day_adj, levels=c(
  "-31","-13","-3","0","0.5","1","1.5","2","3","4","5","6","7","8",
  "11","13","14","16","18","21","23","25","28","35","49","63"))

#take only the days common to all 3 cohorts
good.days = as.factor(c(-3,0,1,2,3,4,5,7,8,11,14,18,21,28))
a1.sub = subset(coefs.meta, Rec_day_adj %in% good.days)

cols=c("goldenrod1","#c46206","#ff7b00","#ffa857","#ffd5ad","#047518","#2a963d","#8bc997",
       "#490073","#9436c9","#c893e6","#02207a","#1c41b0","#6b87db")


#plotting all recovery days through 63
p1 = ggplot(a1.sub, aes(x=Axis.1, y=Axis.2,color=Rec_day_adj,shape=Treatment)) +
  facet_grid(Cohort~Diet) +
  geom_point(size = 3 ,alpha=0.8) + theme_bw() +
  scale_shape_manual(values=c(16,1,17,2)) +
  scale_color_manual(values=cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_text(size = 9, face="bold"), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1)) +
  xlab("PC1 (26.2%)") +ylab("PC2 (12.2%)")
p1
#ggsave(file="FILE.svg", plot=p1, width=7, height=6)



#subset just recovery days through D14
d14 = subset(a1.sub, Rec_day < 16)
p = ggplot(d14, aes(x=Axis.1, y=Axis.2,color=Rec_day_adj,shape=Treatment)) +
  facet_wrap(Diet~.,nrow=2) +
  geom_point(size = 3, alpha=0.9) + theme_bw() +
  scale_shape_manual(values=c(16,1,17,2)) +
  scale_color_manual("Day",values=cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_text(size = 9, face="bold"), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"),
        axis.title=element_blank(),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p
#ggsave(file="FILE.svg", plot=p, width=3.5, height=4.5)



####################################
########### statistics #############
####################################

### PERMANOVAs ###
#results presented in Table S2D

#1. D14 BC dist ~ Abx * Diet
#i.e. at D14, do samples cluster significantly on the basis of ABX, diet, or their interaction?
d14 = subset_samples(ps_gen, Rec_day_adj=="14")
d14.meta = meta(d14)
d14_dist = data.frame(as.matrix(phyloseq::distance(d14, method="bray")))
adonis2(d14_dist ~ d14.meta$Abx*d14.meta$Diet, data = d14.meta, perm=9999)

#posthoc comparisons across treatment groups
posthoc1 = pairwise.adonis(d14_dist, d14.meta$Treatment, p.adjust.m="BH")
posthoc1

#2. D28 BC dist ~ Abx * Diet
d28 = subset_samples(ps_gen, Rec_day_adj=="28")
d28.meta = meta(d28)
d28_dist = data.frame(as.matrix(phyloseq::distance(d28, method="bray")))
adonis2(d28_dist ~ d28.meta$Abx*d28.meta$Diet, data = d28.meta, perm=9999)

#posthoc comparisons across treatment groups
posthoc2 = pairwise.adonis(d28_dist, d28.meta$Treatment, p.adjust.m="BH")
posthoc2


#3. BC dist ~ Diet * Recovery Day
#i.e. for ABX-treated groups, compare across timepoints and treatments
#i.e. do ABX-treated samples return to their pre-ABX baseline across diets?

#subset timepoints to minimize excessive comparisons
keep = c(-3,0,7,14,28)
ps_sub = subset_samples(ps_gen, Rec_day_adj %in% keep)
#take out the pbs controls
ps_sub1 = subset_samples(ps_sub, Abx=="ABX")
ps_sub_meta = meta(ps_sub1)
ps_sub_dist = data.frame(as.matrix(phyloseq::distance(ps_sub1, method="bray")))

adonis2(ps_sub_dist ~ ps_sub_meta$Diet*ps_sub_meta$Rec_day_adj, data = ps_sub_meta, perm=9999)

#dummy variable for day*diet
ps_sub_meta$day_diet = paste0(ps_sub_meta$Diet, "_", ps_sub_meta$Rec_day_adj)

#posthoc comparisons across day/diet combinations
posthoc3 = pairwise.adonis(ps_sub_dist, ps_sub_meta$day_diet, p.adjust.m="BH")
posthoc3






### Beta Diversity ###
#i.e. BC distance from ABX groups on each diet to respective PBS controls at D0, D7, D14, and D28
#results analyzed in Prism and presented in Table S2E


#GENERAL OUTLINE:
#1. Subset a timepoint
#2. For each ABX sample, calculate BC distance to all PBS samples on same diet
#3. Average the BC distance to PBS controls


#collect the genus-level BC distance matrix for all of the samples
ps_gen_meta = meta(ps_gen)
ps_gen_dist = data.frame(as.matrix(phyloseq::distance(ps_gen, method="bray")))


## ABX vs PBS at D0 ##

#subset D0 and each treatment group separately
d0 = subset(ps_gen_meta, Rec_day_adj=="0")
rc.abx = row.names(d0[d0$Treatment=="RC-ABX",])
rc.pbs = row.names(d0[d0$Treatment=="RC-PBS",])
wd.abx= row.names(d0[d0$Treatment=="WD-ABX",])
wd.pbs = row.names(d0[d0$Treatment=="WD-PBS",])

#empty vector to collect results
d0.rc.BC=NULL

#for each RC-ABX sample i at D0...
for(i in 1:length(rc.abx)){
  tmp = NULL #dummy variable to hold results
  
  #for each RC-PBS sample j at D14...
  for(j in 1:length(rc.pbs)){
    
    #calculate and store BC distance from i to all j's
    #i.e. BC distance from RC-ABX[i] to all RC-PBS[j] samples
    tmp = c(tmp,ps_gen_dist[rc.abx[i],rc.pbs[j]])
  }
  
  #Average the BC distances calculated above and store in a new vector
  d0.rc.BC = c(d0.rc.BC,mean(tmp))
}

#same process for WD-ABX vs WD-PBS
d0.wd.BC=NULL
for(i in 1:length(wd.abx)){
  tmp = NULL
  for(j in 1:length(wd.pbs)){
    tmp = c(tmp,ps_gen_dist[wd.abx[i],wd.pbs[j]])
  }
  d0.wd.BC = c(d0.wd.BC,mean(tmp)) 
}


## ABX vs PBS at D7 ##
#same process as above 

d7 = subset(ps_gen_meta, Rec_day_adj=="7")
rc.abx = row.names(d7[d7$Treatment=="RC-ABX",])
rc.pbs = row.names(d7[d7$Treatment=="RC-PBS",])
wd.abx= row.names(d7[d7$Treatment=="WD-ABX",])
wd.pbs = row.names(d7[d7$Treatment=="WD-PBS",])

d7.rc.BC=NULL
for(i in 1:length(rc.abx)){
  tmp = NULL
  for(j in 1:length(rc.pbs)){
    tmp = c(tmp,ps_gen_dist[rc.abx[i],rc.pbs[j]])
  }
  d7.rc.BC = c(d7.rc.BC,mean(tmp))
}

d7.wd.BC=NULL
for(i in 1:length(wd.abx)){
  tmp = NULL
  for(j in 1:length(wd.pbs)){
    tmp = c(tmp,ps_gen_dist[wd.abx[i],wd.pbs[j]])
  }
  d7.wd.BC = c(d7.wd.BC,mean(tmp))
}



## ABX vs PBS at D14 ##
#same process as above
d14 = subset(ps_gen_meta, Rec_day_adj=="14")
rc.abx = row.names(d14[d14$Treatment=="RC-ABX",])
rc.pbs = row.names(d14[d14$Treatment=="RC-PBS",])
wd.abx= row.names(d14[d14$Treatment=="WD-ABX",])
wd.pbs = row.names(d14[d14$Treatment=="WD-PBS",])

d14.rc.BC=NULL

for(i in 1:length(rc.abx)){
  tmp = NULL 
  for(j in 1:length(rc.pbs)){
    tmp = c(tmp,ps_gen_dist[rc.abx[i],rc.pbs[j]])
  }
  d14.rc.BC = c(d14.rc.BC,mean(tmp)) #mean(tmp)
}

#do the same thing for WD-ABX vs WD-PBS
d14.wd.BC=NULL
for(i in 1:length(wd.abx)){
  tmp = NULL
  for(j in 1:length(wd.pbs)){
    tmp = c(tmp,ps_gen_dist[wd.abx[i],wd.pbs[j]])
  }
  d14.wd.BC = c(d14.wd.BC,mean(tmp))
}


## ABX vs PBS at D28 ##
#same process as above
d28 = subset(ps_gen_meta, Rec_day_adj=="28")
rc.abx = row.names(d28[d28$Treatment=="RC-ABX",])
rc.pbs = row.names(d28[d28$Treatment=="RC-PBS",])
wd.abx= row.names(d28[d28$Treatment=="WD-ABX",])
wd.pbs = row.names(d28[d28$Treatment=="WD-PBS",])

d28.rc.BC=NULL
for(i in 1:length(rc.abx)){
  tmp = NULL
  for(j in 1:length(rc.pbs)){
    tmp = c(tmp,ps_gen_dist[rc.abx[i],rc.pbs[j]])
  }
  d28.rc.BC = c(d28.rc.BC,mean(tmp))
}

d28.wd.BC=NULL
for(i in 1:length(wd.abx)){
  tmp = NULL
  for(j in 1:length(wd.pbs)){
    tmp = c(tmp,ps_gen_dist[wd.abx[i],wd.pbs[j]])
  }
  d28.wd.BC = c(d28.wd.BC,mean(tmp))
}


#these vectors were imported into prism for two-way ANOVA and plotting

