library(tidyverse)
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(gcookbook)
library(scales)
library(RColorBrewer)

#this script calculates and plots the RELATIVE ABUNDANCE of sig. enriched KOs 
#across timepoints and diets. 
#plots show the abundance of the KOs at D-3 pre-ABX next to their abundance at
#the later timepoint at which they are significantly enriched

setwd("~/Documents/Chang-Bergelson Lab/HF Diet Abx/Aim 1/Metagenomics/comm_analysis/")

############################################################################
### 1. Data processing ####
### need to take the enriched KOs and grab all their relative abundances
### the output of this part is availableat 02.data_files/meanRelAbun_enr.csv
############################################################################

counts = read.csv("counts_ko_sum.csv") #raw abundances
#ko = read.csv("ko.mod4.csv") #KEGG mapping
#kcat = dplyr::select(ko, KO, kcat.mod) #KEGG cat to KO
#kcat = unique(kcat) #no duplicates
#ksys = dplyr::select(ko, c("ksys.mod","kfam.mod","kcat.mod")) #full hierarchy
#ksys = unique(ksys) #no duplicates

#convert KO abundances to relative abundances
#rel.cov = counts
#for(i in 2:ncol(rel.cov)){
#  tot = sum(rel.cov[,i])
#  rel.cov[,i] = rel.cov[,i]/tot
#}

#clean up dataframe
#row.names(rel.cov) = rel.cov$KO
#rel.cov1 = rel.cov[,-1] #get rid of the column for KO

#convert to matrix
#ko.m = as.matrix(rel.cov1)

#to import to phyloseq to tidy
#meta = read.csv("meta.csv")
#row.names(meta) = meta$SeqID

#OTU = otu_table(ko.m, taxa_are_rows = TRUE)
#META = sample_data(meta)
#ps = phyloseq(OTU,META)

#sample_data(ps)$Rec_day_adj = factor(sample_data(ps)$Rec_day_adj, levels=c("-13","-3","2","4","14","28"))

#subset only the mice for which we have all the timepoints
#complete = c(1507, 1512, 1946, 1953, 5733, 5735)
#ps.cpt = subset_samples(ps, Mouse %in% complete)

#subset RC-ABX and WD-ABX groups
#rc = subset_samples(ps.cpt, Treatment=="RC-ABX")
#wd = subset_samples(ps.cpt, Treatment=="WD-ABX")

#put OTU table into long/tidy form
#rc_df = psmelt(rc)
#colnames(rc_df)[colnames(rc_df)=="OTU"] = "KO"

#import the list of sig. enriched KOs
#d2.DEGs = read.csv("DEG_vs_preABX/rc.d2.DEG.csv")
#d2.enr = subset(d2.DEGs, log2FoldChange > 0) #only take the enriched ones, not depleted

#subset the dataframe for the enriched KOs
#rc.d2 = subset(rc_df, KO %in% d2.enr$KO)

#map to KEGG cat level
#rc.d2.ko = left_join(rc.d2, kcat, multiple="first")

#label as enriched at D2
#rc.d2.ko$day = "D2"

#import/subset enriched KOs for RC D4
#d4.DEGs = read.csv("deseq_DEGs_byDay/rc.d4.DEG.csv")
#d4.enr = subset(d4.DEGs, log2FoldChange > 0)
#rc.d4 = subset(rc_df, KO %in% d4.enr$KO)
#rc.d4.ko = left_join(rc.d4, kcat, multiple="first")
#rc.d4.ko$day = "D4"

#merge D2 and D4 enriched dataframes
#rc.comb = rbind(rc.d2.ko, rc.d4.ko)
#rc.comb$day = factor(rc.comb$day, levels=c("D2","D4"))


#same process for WD
#wd_df = psmelt(wd)
#colnames(wd_df)[colnames(wd_df)=="OTU"] = "KO"

#wd.d2.DEGs = read.csv("DEG_vs_preABX/wd.d2.DEG.csv")
#wd.d2.enr = subset(wd.d2.DEGs, log2FoldChange > 0)
#wd.d2 = subset(wd_df, KO %in% wd.d2.enr$KO)
#wd.d2.ko = left_join(wd.d2, kcat, multiple="first")
#wd.d2.ko$day = "D2"

#wd.d14.DEGs = read.csv("DEG_vs_preABX/wd.d14.DEG.csv")
#wd.d14.enr = subset(wd.d14.DEGs, log2FoldChange > 0)
#wd.d14 = subset(wd_df, KO %in% wd.d14.enr$KO)
#wd.d14.ko = left_join(wd.d14, kcat, multiple="first")
#wd.d14.ko$day = "D14"

#wd.d28.DEGs = read.csv("DEG_vs_preABX/wd.d28.DEG.csv")
#wd.d28.enr = subset(wd.d28.DEGs, log2FoldChange > 0)
#wd.d28 = subset(wd_df, KO %in% wd.d28.enr$KO)
#wd.d28.ko = left_join(wd.d28, kcat, multiple="first")
#wd.d28.ko$day = "D28"

#wd.comb = rbind(wd.d2.ko, wd.d14.ko, wd.d28.ko)
#wd.comb$day = factor(wd.comb$day, levels=c("D2","D14", "D28"))
#wd.comb = left_join(wd.comb, ksys, multiple="first")


#combine the WD and RC datasets
#comb.all = rbind(wd.comb, rc.comb)
#comb.all$Rec_day_adj = factor(comb.all$Rec_day_adj, levels=c("-13","-3","2","4","14","28"))

#average the relative abundances within each diet/timepoint across mice
#mean.comb = comb.all %>% group_by(Treatment, day, Rec_day_adj, KO) %>%
#  summarize(meanRelAbun = mean(Abundance))

#map to KEGG category and system
#mean.comb1 = left_join(mean.comb, kcat, multiple="first")
#mean.comb2 = left_join(mean.comb1, ksys,multiple="first")

#write.csv(mean.comb2, "meanRelAbun_enr.csv",row.names=FALSE)





############################################################################
### 2. Plot the enriched datasets
############################################################################

#this file is generated above, or available in the 02.data_files/ directory
mean.comb2 = read.csv("meanRelAbun_enr.csv")
mean.comb2$kfam.mod <- factor(mean.comb2$kfam.mod,
                                  levels = unique(mean.comb2$kfam.mod), ordered = TRUE)

#take out pre-diet-acclimation timepoint
mean.comb2 = subset(mean.comb2, Rec_day_adj != "-13")

mean.comb2$Rec_day_adj = factor(mean.comb2$Rec_day_adj, levels=c("-3","2","4","14","28"))

cols=c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#ff9900", "grey70",  
       "#1B9E77","#D95F02", "#f78b83","#7570B3",  "#66A61E", "#E6AB02", "#A6761D", "#A6CEE3",
       "#B2DF8A", "#33A02C", "#FB9A99","#FDBF6F" ,"#E31A1C","#FFFF99", "turquoise3", "#FF7F00", 
       "#7FC97F", "#BEAED4", "#FDC086")


rc.comb = subset(mean.comb2, Treatment=="RC-ABX")

#sum the rel abundances for all the KOs that map to one KEGG family
rc.comb1 = rc.comb %>% group_by(day, Rec_day_adj, kfam.mod) %>% 
  summarize(kfam.relAbun = sum(meanRelAbun))


#### just the D2 enriched ones
rc.2 = subset(rc.comb1, day=="D2") %>%
  subset(Rec_day_adj=="-3" | Rec_day_adj=="2")

#redo the colors again
rc.cols2 = c("#7FC97F","#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#ff9900", "grey70",  
            "#D95F02", "#7570B3", "#E6AB02", "#A6761D", "#A6CEE3",
            "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "turquoise3", "#FF7F00", 
            "#7FC97F", "#BEAED4", "#FDC086")


#Plots the rel abundance of KOs enriched at D2 relative to pre-ABX D-3 for RC
#aka Figure S2J
p = ggplot(data=rc.2, aes(x=Rec_day_adj, y=kfam.relAbun, fill=kfam.mod)) + 
  geom_bar(stat="identity",position="stack") +
  scale_fill_manual(values=rc.cols2) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5)),
         fill=guide_legend(override.aes = list(size = 0.5), ncol=1)) +
  facet_grid(.~day) +
  scale_y_continuous(limits=c(0,0.23),breaks=c(seq(0,0.2,.05))) +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_line(color="grey90"),
    panel.grid.minor = element_line(color="grey90"),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12,color="black"),
        #axis.ticks.y = element_blank(),
        #strip.text.x = element_text(face="bold",size=10),
        strip.text.x = element_blank(),
        #strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        #legend.title = element_text(face='bold'),
        #legend.text  = element_text(size = 7),
        #legend.key.size = unit(0.6, "lines")) +
        legend.position="none") +
  guides(fill=guide_legend(title="KEGG Family",ncol=1))
p
#ggsave(file="RC_enr_kfam_d2.svg", plot=p, width=1.8, height=3)


### just the D4 enriched ones
rc.4 = subset(rc.comb1, day=="D4") %>%
  subset(Rec_day_adj=="-3" | Rec_day_adj=="4")

#redo the colors again
rc.cols4 = c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#ff9900", "grey70",  
            "#1B9E77","#D95F02", "#f78b83","#7570B3",
            "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C","#FFFF99", "turquoise3", "#FF7F00", 
            "#7FC97F", "#BEAED4", "#FDC086")


#Plots the rel abundance of KOs enriched at D4 relative to pre-ABX D-3 for RC
#aka Figure S2K
p1 = ggplot(data=rc.4, aes(x=Rec_day_adj, y=kfam.relAbun, fill=kfam.mod)) + 
  geom_bar(stat="identity",position="stack") +
  scale_fill_manual(values=rc.cols4) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5)),
         fill=guide_legend(override.aes = list(size = 0.5), ncol=1)) +
  facet_grid(.~day) +
  scale_y_continuous(limits=c(0,0.23),breaks=c(seq(0,0.2,.05))) +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_line(color="grey90"),
    panel.grid.minor = element_line(color="grey90"),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12,color="black"),
        #axis.ticks.y = element_blank(),
        #strip.text.x = element_text(face="bold",size=10),
        strip.text.x = element_blank(),
        #strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        #legend.title = element_text(face='bold'),
        #legend.text  = element_text(size = 7),
        #legend.key.size = unit(0.6, "lines")) +
        legend.position="none") +
  guides(fill=guide_legend(title="KEGG Family",ncol=1))
p1
#ggsave(file="RC_enr_kfam_d4.svg", plot=p1, width=1.8, height=3)



##same for WD

wd.comb = subset(mean.comb2, Treatment=="WD-ABX")
wd.comb1 = wd.comb %>% group_by(day, Rec_day_adj, kfam.mod) %>% 
  summarize(kfam.relAbun = sum(meanRelAbun))
wd.comb1$day = factor(wd.comb1$day, levels=c("D2","D14","D28"))

#### just the D2 enriched ones

wd.2 = subset(wd.comb1, day=="D2") %>%
  subset(Rec_day_adj=="-3" | Rec_day_adj=="2")

#redo the colors again
wd.cols2=c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#ff9900", "grey70",  
          "#1B9E77","#D95F02", "#7570B3",  "#66A61E", "#A6761D",
          "#B2DF8A", "#33A02C", "#FB9A99","#FDBF6F" ,"#E31A1C", "turquoise3", "#FF7F00", 
          "#7FC97F", "#BEAED4", "#FDC086")


#Plots the rel abundance of KOs enriched at D2 relative to pre-ABX D-3 for WD
#aka Figure S2M
p3 = ggplot(data=wd.2, aes(x=Rec_day_adj, y=kfam.relAbun, fill=kfam.mod)) + 
  geom_bar(stat="identity",position="stack") +
  scale_fill_manual(values=wd.cols2) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5)),
         fill=guide_legend(override.aes = list(size = 0.5), ncol=1)) +
  facet_grid(.~day) +
  scale_y_continuous(limits=c(0,0.23),breaks=c(seq(0,0.2,.05))) +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_line(color="grey90"),
    panel.grid.minor = element_line(color="grey90"),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12,color="black"),
        #axis.ticks.y = element_blank(),
        #strip.text.x = element_text(face="bold",size=10),
        strip.text.x = element_blank(),
        #strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        #legend.title = element_text(face='bold'),
        #legend.text  = element_text(size = 7),
        #legend.key.size = unit(0.6, "lines")) +
        legend.position="none") +
  guides(fill=guide_legend(title="KEGG Family",ncol=1))
p3
#ggsave(file="WD_enr_kfam_d2.svg", plot=p3, width=1.8, height=3)



### just the D4 enriched ones
wd.14 = subset(wd.comb1, day=="D14") %>%
  subset(Rec_day_adj=="-3" | Rec_day_adj=="14")

#redo the colors again
wd.cols14=c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#ff9900", "grey70",  
          "#1B9E77","#D95F02", "#f78b83","#7570B3",  "#66A61E", "#A6761D",
          "#B2DF8A", "#33A02C", "#FB9A99","#FDBF6F" ,"#E31A1C","#FFFF99", "#FF7F00", 
          "#7FC97F", "#BEAED4", "#FDC086")

#Plots the rel abundance of KOs enriched at D14 relative to pre-ABX D-3 for WD
#aka Figure S2N
p4 = ggplot(data=wd.14, aes(x=Rec_day_adj, y=kfam.relAbun, fill=kfam.mod)) + 
  geom_bar(stat="identity",position="stack") +
  scale_fill_manual(values=wd.cols14) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5)),
         fill=guide_legend(override.aes = list(size = 0.5), ncol=1)) +
  facet_grid(.~day) +
  scale_y_continuous(limits=c(0,0.23),breaks=c(seq(0,0.2,.05))) +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_line(color="grey90"),
    panel.grid.minor = element_line(color="grey90"),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12,color="black"),
        #axis.ticks.y = element_blank(),
        #strip.text.x = element_text(face="bold",size=10),
        strip.text.x = element_blank(),
        #strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        #legend.title = element_text(face='bold'),
        #legend.text  = element_text(size = 7),
        #legend.key.size = unit(0.6, "lines")) +
        legend.position="none") +
  guides(fill=guide_legend(title="KEGG Family",ncol=1))
p4
#ggsave(file="WD_enr_kfam_d4.svg", plot=p4, width=1.8, height=3)
