library(tidyverse)
library(phyloseq)
library(microbiome)
library(RColorBrewer)

setwd("/DIRECTORY/")

meta = read.csv("meta.csv")
row.names(meta) = meta$SeqID

#file to match up gene caller IDs with KOs
kofams = read_table("wd_kofams.txt")
kofams = dplyr::select(kofams, c("gene_callers_id","accession")) 
colnames(kofams)[2] = "KO"

kofams$gene_callers_id = as.character(kofams$gene_callers_id)

counts <-read.csv("wd_all_cov.csv",row.names=1)
counts1 = counts[colnames(counts) %in% meta$SeqID]

counts1$gene_callers_id = row.names(counts1)

#put the gene call counts in long form
counts.tidy = gather(counts1, key="SeqID",value="Count", -gene_callers_id)

#match up gene calls with KOs
counts.tidy.ko = left_join(counts.tidy, kofams, multiple="all", relationship="many-to-many")

#filter to include only genes that are present
counts.ko.non0 = subset(counts.tidy.ko, Count > 0)

#filter out the genes that don't have a known KO
counts.ko.non0 = counts.ko.non0[!is.na(counts.ko.non0$KO),]

#merge with metadata
counts.ko.meta = left_join(counts.ko.non0, meta)

genesPerKO = counts.ko.meta %>% group_by(KO, Treatment, Rec_day_adj, Mouse) %>% 
  summarize(nGenes = length(unique(gene_callers_id)))

genesPerKO$Treatment = factor(genesPerKO$Treatment, levels=c("WD-ABX","RC-ABX"))

gpk.wide = spread(genesPerKO, key=KO, value=nGenes)
gpk.wide[is.na(gpk.wide)] = 0
genesPerKO = gather(gpk.wide, key=KO, value=nGenes, -Treatment,-Rec_day_adj,-Mouse)

#for each mouse, calculate its mean genes per KO for each timepoint
#this is a measure of functional redundancy
mean_genesPerKO = genesPerKO %>%
  group_by(Treatment, Rec_day_adj, Mouse) %>%
  summarize(meanPerKO = mean(nGenes))

#average functional redundancy across mice at each timepoint
fr.summ = mean_genesPerKO %>% group_by(Rec_day_adj, Treatment) %>% 
  summarise(
    mean_fr = mean(meanPerKO),
    sd = sd(meanPerKO, na.rm=TRUE))

#write.csv(mean_genesPerKO, "func_redund_raw.csv",row.names=FALSE)
#write.csv(fr.summ, "func_redund_mean.csv",row.names=FALSE)

# ^^^ these files are used in 04.func_richness_redund.pzfx



#plot to see if the KOs with high functional redundancy at D28 are just
#the ones that started with a lot of functional redundancy
#i.e. compare initial (D-3) vs final (D28) FR
#This produces Figure S2B, Table S3E

fr.wide = spread(genesPerKO, key=Rec_day_adj, value=nGenes)
fr.wide[is.na(fr.wide)] = 0

fr.rc = subset(fr.wide, Treatment=="RC-ABX")
fr.wd = subset(fr.wide, Treatment=="WD-ABX")

fr.wide$Treatment = factor(fr.wide$Treatment, levels=c("RC-ABX","WD-ABX"))

#Figure S2B:
p = ggplot(fr.wide, aes(x=`-3`,y=`28`, color=Treatment)) + geom_point(size=2,alpha=0.7) + 
  xlab("Gene Calls per KO: D-3") + 
  ylab("Gene Calls per KO: D28") +
  scale_y_log10() + scale_x_log10() +
  #scale_y_continuous(breaks=c(0,500,1000),limits=c(0,1200)) +
  scale_color_manual(values=c("cornflowerblue","salmon2")) +
  geom_abline(intercept=0,slope=1,lty=2,colour="darkgray") + theme_bw() +
  theme(#panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_text(size = 9, face="bold"), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1)) +
  facet_wrap(Treatment ~.)
p
#ggsave(file=initial_final_FR.png, plot=p, width=7, height=4)

#subset by diet and do a linear regression

fr.rc = subset(fr.wide, Treatment=="RC-ABX")
fr.wd = subset(fr.wide, Treatment=="WD-ABX")

#only analyze the ones that came back at all
fr.rc = subset(fr.rc, `28` > 0)
fr.wd = subset(fr.wd, `28` > 0)


rc.lm = lm(fr.rc$`28`~fr.rc$`-3`)
summary(rc.lm)

wd.lm = lm(fr.wd$`28`~fr.wd$`-3`)
summary(wd.lm)




##################################################################################
########### let's look at which KOs have the most/least recoverable FR ##########

#average mean genes per KO across mice within a treatment group
#i.e., genesPerKO is the average genes per KO for all the gene calls within one mouse
#genes by ko is the average of the mean genes per KO across all the mice within a treatment/timepoint
genes_by_ko = genesPerKO %>% group_by(KO, Rec_day_adj, Treatment) %>% 
  summarize(meanGPK = mean(nGenes))


#gather the pre-abx timepoint so we can calculate FR as percentage of initial value at D-3
d3 = subset(genes_by_ko, Rec_day_adj==-3)
d3$ko.trt = paste(d3$KO, d3$Treatment, sep="_")
colnames(d3)[4] = "pre"
d3 = dplyr::select(d3, c(pre, ko.trt))
new.d3 = data.frame(pre = d3$pre, ko.trt = d3$ko.trt)

#dummy variable for merging with pre-abx timepoint
genes_by_ko$ko.trt = paste(genes_by_ko$KO, genes_by_ko$Treatment, sep="_")
#merge with pre-abx timepoint
genes_by_ko = left_join(genes_by_ko, new.d3)
#calculate FR as percentage of initial genes per KO
genes_by_ko$perc_rec = genes_by_ko$meanGPK/genes_by_ko$pre

d14 = subset(genes_by_ko, Rec_day_adj==14)
d14$perc_rec[d14$perc_rec==Inf] = NA
d14$perc_rec[is.nan(d14$perc_rec)] = NA
d14.rc = subset(d14, Treatment=="RC-ABX")
d14.wd = subset(d14, Treatment=="WD-ABX")

#histogram showing how often genes recovered a given percentage of FR
hist(d14.rc$perc_rec, breaks=40,xlim=range(0,2),xlab="% FR Recovery",ylab="nKOs")
hist(d14.wd$perc_rec, breaks=400,xlim=range(0,2),xlab="% FR Recovery",ylab="nKOs")

#this gathers the set of KOs that recovered well (>75% recovery) or poorly (<25%)
top25 = subset(d14, perc_rec > 0.75)
bottom25 = subset(d14, perc_rec < .25)





#now map the good and poor recoverers by KEGG system
#This produces Figure 2B

ko = read.csv("ko.mod4.csv")
ko.neat = dplyr::select(ko, c(KO, kcat.mod, kfam.mod, ksys.mod))


top25.ko = left_join(top25, ko.neat, multiple="all", relationship="many-to-many")
bottom25.ko = left_join(bottom25, ko.neat, multiple="all", relationship="many-to-many")


#tally up how many of the genes that recovered well/poorly mapped to each category
t25.tally = top25.ko %>% group_by(Treatment, ksys.mod) %>% tally()
b25.tally = bottom25.ko %>% group_by(Treatment, ksys.mod) %>% tally()

t25.tally$top = "Top25"
b25.tally$top = "Bottom25"

all.tally = rbind(t25.tally, b25.tally)
all.tally = all.tally[!is.na(all.tally$ksys.mod),]


qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cols=col_vector[1:24]
cols[12] = "#f78b83"
cols[7] = "#ff9900"
cols[8] = "grey70"
cols=rep(cols, 20)

p = ggplot(data = all.tally, aes(x=top, fill=ksys.mod, weight=n)) + 
  geom_bar(position="stack") + #theme(legend.position="none") +
  scale_fill_manual(values=cols) +ylab("nKOs") +
  facet_wrap(.~Treatment) + theme_bw() +
  theme(#panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_text(size = 9, face="bold"), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"),
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p
#ggsave(file="FILE.svg", plot=p, width=4.5, height=3)

