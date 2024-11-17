library(tidyverse)
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(gcookbook)
library(scales)
library(RColorBrewer)

setwd("/DIRECTORY/")


####  start with GENES LEVEL, then go to codes   ###
otu = read.csv("wd_merged_otu.csv", row.names = 1)
meta = read.csv("wd_meta_merged.csv", row.names=1)
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

OTU = otu_table(otu_matrix, taxa_are_rows = TRUE)
META = sample_data(meta)
TAX = tax_table(tax_matrix)
TREE = read_tree("wd_merged_tree.nwk")
ps = phyloseq(OTU,META,TAX, TREE)

a2 = subset_samples(ps, Aim=="A2")

sample_data(a2)$Rec_day_adj = factor(sample_data(a2)$Rec_day_adj, levels=c(
  "-31","-3","0","7","14","28"))

#filter out less than 5000 reads/sample
ps1 = prune_samples(sample_sums(a2)>=5000, a2)
#transform to rel abundance
ps2 = transform_sample_counts(ps1, function(x) x / sum(x) )

#subset D14
d14 = subset_samples(ps2, Rec_day_adj=="14")

#put it in tidy/long form
d14_df = psmelt(d14)

#sum rel abundances at family level
d14_fam = d14_df %>% group_by(Sample, Kingdom, Phylum, Class, Order, Family) %>% 
  summarise(Abun = sum(Abundance))

meta$Sample = rownames(meta)

#merge tidy otu abundances with metadata
d14_fam1 = left_join(d14_fam, meta, by="Sample")

#average rel abundances at family level across mice for each timepoint/treatment group
d14_mean = d14_fam1 %>% group_by(Treatment, Rec_day_adj, Kingdom, Phylum, Class, Order, Family) %>%
  summarise(meanAbun = mean(Abun))

#rename phyla with < 1% abundance
d14_mean$Family[d14_mean$meanAbun < 0.01] <- "< 1% abund."
d14_mean$Phylum[d14_mean$meanAbun < 0.01] <- "< 1% abund."

#rename the NAs to "unknown"
d14_mean$Family[d14_mean$Family==""] = "Unknown"
d14_mean$Phylum[d14_mean$Phylum==""] = "Unknown"

adj_cols=rep(NA, 18)
adj_cols[1] = "orchid" #<1% abund
adj_cols[2] = "green4" #alcaligenaceae
adj_cols[3] = "cyan3" #Bacteroidaceae
adj_cols[4] = "#BF5B17" #Bacteroidales S24-7
adj_cols[5] = "#FFFF99" #Bifidobacteriaceae
adj_cols[6] = "#1B9E77" #Clostridiaceae 1
#adj_cols[7] =  "grey"#Corio
adj_cols[7] = "darkgreen" #Enterococcaceae
adj_cols[8] = "orange" #Erysipelo
adj_cols[9] = "salmon" #Lachno
adj_cols[10] = "royalblue4" #Lacto
adj_cols[11] = "#A6CEE3" #Peptostrep
adj_cols[12] = "aquamarine4" #porphyromonad
adj_cols[13] = "darkred" #rikenellaceae
adj_cols[14] = "#a0be6c" #Ruminococcus
adj_cols[15] = "lemonchiffon3" #Streptococc
adj_cols[16] = "darkorange2" #verrucomicrobiaceae
adj_cols[17] = "#6A3D9A" #Xanthomonad
adj_cols[18] = "#CAB2D6" #Unknown


d14_mean$Treatment = factor(d14_mean$Treatment, levels=c(
  "RC-PBS-PBS","WD-PBS-PBS","RC-RC","RC-MICR","RC-ABX-PBS","WD-DIET","WD-RC",
  "WD-WD","WD-MICR","WD-ABX-PBS","RC-DIET","RC-WD"))

#merge all the <1% guys into one category
clean = d14_mean %>% group_by(Treatment, Family) %>% 
  summarize(meanAbun = sum(meanAbun))


p1 = ggplot(data=clean, aes(x=Treatment, y=meanAbun, fill=Family)) + 
  geom_bar(stat="identity",position="stack") +
  scale_fill_manual("Family", values=adj_cols) +
  scale_y_continuous(breaks=c(seq(0,1,0.1)),expand=expansion(mult=c(0.01,0.01))) +
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
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        legend.title = element_text(face='bold'),
        legend.position = "none",
        strip.text.x = element_text(face="bold",size=12),
        strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p1
#ggsave(file="FILE.svg", plot=p1, width=5, height=2.5)




#now D28
d28 = subset_samples(ps2, Rec_day_adj=="28")

#put it in tidy/long form
d28_df = psmelt(d28)

#sum rel abundances at family level
d28_fam = d28_df %>% group_by(Sample, Kingdom, Phylum, Class, Order, Family) %>% 
  summarise(Abun = sum(Abundance))

meta$Sample = rownames(meta)

#merge tidy otu abundances with metadata
d28_fam1 = left_join(d28_fam, meta, by="Sample")

#average rel abundances at family level across mice for each timepoint/treatment group
d28_mean = d28_fam1 %>% group_by(Treatment, Rec_day_adj, Kingdom, Phylum, Class, Order, Family) %>%
  summarise(meanAbun = mean(Abun))

#rename phyla with < 1% abundance
d28_mean$Family[d28_mean$meanAbun < 0.01] <- "< 1% abund."
d28_mean$Phylum[d28_mean$meanAbun < 0.01] <- "< 1% abund."

#rename the NAs to "unknown"
d28_mean$Family[d28_mean$Family==""] = "Unknown"
d28_mean$Phylum[d28_mean$Phylum==""] = "Unknown"

#set color scheme to match other plots
adj_cols = rep(NA, 18)
adj_cols[1] = "orchid" #<1% abund
adj_cols[2] = "green4" #alcaligenaceae
adj_cols[3] = "cyan3" #Bacteroidaceae
adj_cols[4] = "#BF5B17" #Bacteroidales S24-7
adj_cols[5] = "#FFFF99" #Bifidobacteriaceae
adj_cols[6] = "#1B9E77" #Clostridiaceae 1
adj_cols[7] =  "grey"#Corio
adj_cols[8] = "darkgreen" #Enterococcaceae
adj_cols[9] = "orange" #Erysipelo
adj_cols[10] = "salmon" #Lachno
adj_cols[11] = "royalblue4" #Lacto
adj_cols[12] = "#A6CEE3" #Peptostrep
#adj_cols[13] = "aquamarine4" #porphyromonad
adj_cols[13] = "grey40" #rikenellaceae
adj_cols[14] = "#a0be6c" #Ruminococcus
adj_cols[15] = "lemonchiffon3" #Streptococc
adj_cols[16] = "darkorange2" #verrucomicrobiaceae
adj_cols[17] = "#6A3D9A" #Xanthomonad
adj_cols[18] = "#CAB2D6" #Unknown


d28_mean$Treatment = factor(d28_mean$Treatment, levels=c(
  "RC-PBS-PBS","WD-PBS-PBS","RC-RC","RC-MICR","RC-ABX-PBS","WD-DIET","WD-RC",
  "WD-WD","WD-MICR","WD-ABX-PBS","RC-DIET","RC-WD"))

#merge all the <1% guys into one category
clean = d28_mean %>% group_by(Treatment, Family) %>% 
  summarize(meanAbun = sum(meanAbun))


p2 = ggplot(data=clean, aes(x=Treatment, y=meanAbun, fill=Family)) + 
  geom_bar(stat="identity",position="stack") +
  scale_fill_manual("Family", values=adj_cols) +
  scale_y_continuous(breaks=c(seq(0,1,0.1)),expand=expansion(mult=c(0.01,0.01))) +
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
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        legend.title = element_text(face='bold'),
        #legend.position = "none",
        strip.text.x = element_text(face="bold",size=12),
        strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p2
#ggsave(file="FILE.svg", plot=p2, width=5, height=2.5)



