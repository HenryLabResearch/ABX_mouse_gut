library(tidyverse)
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(gcookbook)
library(scales)
library(RColorBrewer)

setwd("/DIRECTORY/")

#import data
otu = read.csv("wd_merged_otu.csv", header=TRUE, row.names = 1)
sample_data = read.csv("wd_meta_merged.csv", row.names=1)
tax = read.csv("wd_taxonomy.csv", row.names = 1, header=TRUE)


#Filter out mitochondrial seqs
mito = row.names(tax)[tax$Family=="Mitochondria"]
otu.noMito = otu[!row.names(otu) %in% mito,]
tax.noMito = tax[!row.names(tax) %in% mito,]

#convert to matrix
otu_matrix = as.matrix(otu.noMito)
tax_matrix = as.matrix(tax.noMito)

#convert to phyloseq objects
OTU = otu_table(otu_matrix, taxa_are_rows = TRUE)
META = sample_data(sample_data)
TAX = tax_table(tax_matrix)
ps = phyloseq(OTU,META,TAX)
sample_data(ps)$Rec_day_adj = as.factor(sample_data(ps)$Rec_day_adj)
sample_data(ps)$Rec_day_adj = factor(sample_data(ps)$Rec_day_adj, levels=c(
  "-31","-13","-3","0","0.5","1","1.5","2","3","4","5","6",
  "7","8","11","13","14","16","18","21","23","25","28","35","49","63"))

#filter out less than 5000 reads/sample
ps1 = prune_samples(sample_sums(ps)>=5000, ps)
#transform to rel abundance
ps2 = transform_sample_counts(ps1, function(x) x / sum(x) )

#subset a1 (initial WD-resilience experiements)
a1 = subset_samples(ps2, Aim=="A1")
a1.c1 = subset_samples(a1, Cohort=="C1")
a1.c2 = subset_samples(a1, Cohort=="C2")
a1.c3 = subset_samples(a1, Cohort=="C3")

#just C1 first
c1_df = psmelt(a1.c1)

#sum abundances at the family level
c1_fam = c1_df %>% group_by(Sample, Kingdom, Phylum, Class, Order, Family) %>% 
  summarise(Abun = sum(Abundance))

sample_data$Sample = rownames(sample_data)

#join family level abundances w metadata
c1_fam1 = left_join(c1_fam, sample_data, by="Sample")

#average across mice
c1_mean = c1_fam1 %>% group_by(Treatment, Rec_day_adj, Kingdom, Phylum, Class, Order, Family) %>%
  summarise(meanAbun = mean(Abun))

#rename phyla with < 1% abundance
c1_mean$Family[c1_mean$meanAbun < 0.01] <- "< 1% abund."
c1_mean$Phylum[c1_mean$meanAbun < 0.01] <- "< 1% abund."

#rename the NAs to "unknown"
c1_mean$Family[c1_mean$Family==""] = "Unknown"
c1_mean$Phylum[c1_mean$Phylum==""] = "Unknown"

#make numeric recovery days into a factor
c1_mean$Rec_day_adj = factor(c1_mean$Rec_day_adj, levels=c(
  "-31","-13","-3","0","0.5","1","1.5","2","3","4","5","6",
  "7","8","11","13","14","16","18","21","23","25","28","35","49","63"))

#pick contrasting colors
cols[1] = "orchid" #<1% abund
cols[2] = "green4" #alcaligenaceae
cols[3] = "#BEAED4" #anaeroplasma
cols[4] = "cyan3" #Bacteroidaceae
cols[5] = "#BF5B17" #Bacteroidales S24-7
cols[6] = "#FFFF99" #Bifidobacteriaceae
cols[7] = "#386CB0" #Carnobacteriaceae
cols[8] = "#1B9E77" #Clostridiaceae 1
cols[9] =  "#FDC086"#Clostridiales vadin
cols[10] = "darkred" #Desulfo
#cols[11] = "deeppink" #Enterobacteriaceae
cols[11] = "darkgreen" #Enterococcaceae
cols[12] = "orange" #Erysipelo
cols[13] = "salmon" #Lachno
cols[14] = "royalblue4" #Lacto
cols[15] = "red2" #Morax
cols[16] = "#A6CEE3" #Peptostrep
cols[17] = "#1F78B4" #Planococcaceae
cols[18] = "aquamarine4" #porphyromonad
cols[19] = "goldenrod" #Prevotellaceae
cols[20] = "#666666" #Rikenellaceae
cols[21] = "#a0be6c" #Ruminococcus
cols[22] = "#ffe600" #Staphylococc
cols[23] = "lemonchiffon3" #Streptococc
cols[24] = "#6A3D9A" #Xanthomonad
cols[25] = "#CAB2D6" #Unknown

#we originally had D0.5 and D1.5 but cut them out for visual clarity
no.5 = subset(c1_mean, Rec_day_adj != "0.5" & Rec_day_adj != "1.5") %>%
  subset(Rec_day_adj != "-13")

#plot C1
p1 = ggplot(data=no.5, aes(x=Rec_day_adj, y=meanAbun, fill=Family)) + 
  geom_bar(stat="identity",position="stack") +
  facet_wrap(Treatment~., nrow=2) +
  scale_fill_manual("Family", values=cols) +
  scale_y_continuous(breaks=c(seq(0,1,0.1))) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5)),
         fill=guide_legend(override.aes = list(size = 0.5), ncol=1)) +
  theme(legend.title = element_text(size = 7), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.4, "lines")) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        legend.title = element_text(face='bold'),
        strip.text.x = element_text(face="bold",size=12),
        strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p1
#ggsave(file="FILE.svg", plot=p1, width=8, height=7)




#now cohort 2
c2_df = psmelt(a1.c2)

#sum abundances by family
c2_fam = c2_df %>% group_by(Sample, Kingdom, Phylum, Class, Order, Family) %>% 
  summarise(Abun = sum(Abundance))

sample_data$Sample = rownames(sample_data)

#merge family-level abundance with metadata
c2_fam1 = left_join(c2_fam, sample_data, by="Sample")

#average family-level abundance across mice
c2_mean = c2_fam1 %>% group_by(Treatment, Rec_day_adj, Kingdom, Phylum, Class, Order, Family) %>%
  summarise(meanAbun = mean(Abun))

#rename phyla with < 1% abundance
c2_mean$Family[c2_mean$meanAbun < 0.01] <- "< 1% abund."
c2_mean$Phylum[c2_mean$meanAbun < 0.01] <- "< 1% abund."

#rename the NAs to "unknown"
c2_mean$Family[c2_mean$Family==""] = "Unknown"
c2_mean$Phylum[c2_mean$Phylum==""] = "Unknown"

#make numeric recovery days into a factor
c2_mean$Rec_day_adj = factor(c2_mean$Rec_day_adj, levels=c(
  "-31","-13","-3","0","0.5","1","1.5","2","3","4","5","6",
  "7","8","11","13","14","16","18","21","23","25","28","35","49","63"))


#color scheme to match the C1 plot
cols[1] = "orchid" #<1% abund
cols[2] = "green4" #alcaligenaceae
cols[3] = "cyan3" #Bacteroidaceae
cols[4] = "#BF5B17" #Bacteroidales S24-7
cols[5] = "#FFFF99" #Bifidobacteriaceae
cols[6] = "#386CB0" #Carnobacteriaceae
cols[7] = "#1B9E77" #Clostridiaceae 1
cols[8] =  "#FDC086"#Clostridiales vadin
cols[9] =  "grey"#Corio
cols[10] = "darkred" #Desulfo
cols[11] = "darkgreen" #Enterococcaceae
cols[12] = "orange" #Erysipelo
cols[13] = "salmon" #Lachno
cols[14] = "royalblue4" #Lacto
cols[15] = "red2" #Morax
cols[16] = "#A6CEE3" #Peptostrep
cols[17] = "#1F78B4" #Planococcaceae
cols[18] = "aquamarine4" #porphyromonad
cols[19] = "purple" #rhizobiaceae
cols[20] = "#666666" #Rikenellaceae
cols[21] = "#a0be6c" #Ruminococcus
cols[22] = "#ffe600" #Staphylococc
cols[23] = "lemonchiffon3" #Streptococc
cols[24] = "#CAB2D6" #Unknown
cols[25] = "darkorange2" #verrucomicrobiaceae


p2 = ggplot(data=c2_mean, aes(x=Rec_day_adj, y=meanAbun, fill=Family)) + 
  geom_bar(stat="identity",position="stack") +
  facet_wrap(Treatment~., nrow=2) +
  scale_fill_manual("Family", values=cols) +
  scale_y_continuous(breaks=c(seq(0,1,0.1))) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5)),
         fill=guide_legend(override.aes = list(size = 0.5), ncol=1)) +
  theme(legend.title = element_text(size = 7), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.4, "lines")) +
  ylab("Mean Relative Abundance") + xlab("Day") +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        legend.title = element_text(face='bold'),
        strip.text.x = element_text(face="bold",size=12),
        strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p2
#ggsave(file="FILE.svg", plot=p2, width=10, height=7)




#now cohort 3
c3_df = psmelt(a1.c3)

#sum abundances at the family level
c3_fam = c3_df %>% group_by(Sample, Kingdom, Phylum, Class, Order, Family) %>% 
  summarise(Abun = sum(Abundance))

sample_data$Sample = rownames(sample_data)

#merge family-level abundances with metadata
c3_fam1 = left_join(c3_fam, sample_data, by="Sample")

#average family-level abundances across mice
c3_mean = c3_fam1 %>% group_by(Treatment, Rec_day_adj, Kingdom, Phylum, Class, Order, Family) %>%
  summarise(meanAbun = mean(Abun))

#rename phyla with < 1% abundance
c3_mean$Family[c3_mean$meanAbun < 0.01] <- "< 1% abund."
c3_mean$Phylum[c3_mean$meanAbun < 0.01] <- "< 1% abund."

#rename the NAs to "unknown"
c3_mean$Family[c3_mean$Family==""] = "Unknown"
c3_mean$Phylum[c3_mean$Phylum==""] = "Unknown"

#make numeric recovery days into a factor
c3_mean$Rec_day_adj = factor(c3_mean$Rec_day_adj, levels=c(
  "-31","-13","-3","0","0.5","1","1.5","2","3","4","5","6",
  "7","8","11","13","14","16","18","21","23","25","28","35","49","63"))

#matching color scheme
cols[1] = "orchid" #<1% abund
cols[2] = "green4" #alcaligenaceae
cols[3] = "cyan3" #Bacteroidaceae
cols[4] = "#BF5B17" #Bacteroidales S24-7
cols[5] = "#FFFF99" #Bifidobacteriaceae
cols[6] = "#1B9E77" #Clostridiaceae 1
cols[7] =  "#FDC086"#Clostridiales vadin
cols[8] =  "grey"#Corio
cols[9] = "darkred" #Desulfo
cols[10] = "deeppink" #Enterobacteriaceae
cols[11] = "darkgreen" #Enterococcaceae
cols[12] = "orange" #Erysipelo
cols[13] = "salmon" #Lachno
cols[14] = "royalblue4" #Lacto
cols[15] = "red2" #Morax
cols[16] = "#A6CEE3" #Peptostrep
cols[17] = "aquamarine4" #porphyromonad
cols[18] = "purple" #rhizobiaceae
cols[19] = "#666666" #Rikenellaceae
cols[20] = "#a0be6c" #Ruminococcus
cols[21] = "#ffe600" #Staphylococc
cols[22] = "lemonchiffon3" #Streptococc
cols[23] = "darkorange2" #verrucomicrobiaceae


p3 = ggplot(data=c3_mean, aes(x=Rec_day_adj, y=meanAbun, fill=Family)) + 
  geom_bar(stat="identity",position="stack") +
  facet_wrap(Treatment~., nrow=2) +
  scale_fill_manual("Family", values=cols) +
  scale_y_continuous(breaks=c(seq(0,1,0.1))) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5)),
         fill=guide_legend(override.aes = list(size = 0.5), ncol=1)) +
  theme(legend.title = element_text(size = 7), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.4, "lines")) +
  ylab("Mean Relative Abundance") + xlab("Day") +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        legend.title = element_text(face='bold'),
        strip.text.x = element_text(face="bold",size=12),
        strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p3
#ggsave(file="FILE.svg", plot=p3, width=10, height=7)















