library(fgsea)
library(tidyverse)
library(svglite)


### This script runs a series of GSEAs at the KEGG Category level to detect 
#enrichment of broader categories of KOs rather than individual differential cov

setwd("/DIRECTORY/")

ko = read.csv("ko.mod4.csv")
ko.list = dplyr::select(ko, c(kcat.mod, KO))

ko.levels = dplyr::select(ko, c(ksys.mod, kfam.mod, kcat.mod))
ko.levels = unique(ko.levels)

#limit the ko.list to ones that are present in at least 1 sample
ko_sum = read.csv("counts_ko_sum.csv")

kcat.lists = split(ko.list$KO,ko.list$kcat.mod) #244 kcats at hand

kcat_gsea = function(filename){
  data = read.csv(filename)
  data.stats = data$stat
  names(data.stats) = data$KO
  
  gsea.Res <- fgseaMultilevel(pathways=kcat.lists, stats=data.stats)
  gsea.Res.sig = subset(gsea.Res, padj< 0.05)
  
  return(list(gsea.Res, gsea.Res.sig))
  
}


rc.d2 = kcat_gsea("DEG_vs_preABX/rc.d2.csv")
rc.d2.all = rc.d2[[1]]
rc.d2_sig = rc.d2[[2]] #16 pathways
rc.d2_sig = rc.d2_sig %>% rename(pathway = "kcat.mod")
rc.d2_sig = left_join(rc.d2_sig, ko.levels)
rc.d2_sig$leadingEdge = as.character(rc.d2_sig$leadingEdge)
#write.csv(rc.d2_sig, "rc.d2.met.gsea.csv", row.names=FALSE)

rc.d4 = kcat_gsea("DEG_vs_preABX/rc.d4.csv")
rc.d4_sig = rc.d4[[2]] #13 pathways
rc.d4_sig = rc.d4_sig %>% rename(pathway = "kcat.mod")
rc.d4_sig = left_join(rc.d4_sig, ko.levels)
rc.d4_sig$leadingEdge = as.character(rc.d4_sig$leadingEdge)
#write.csv(rc.d4_sig, "rc.d4.met.gsea.csv", row.names=FALSE)

rc.d14 = kcat_gsea("DEG_vs_preABX/rc.d14.csv")
rc.d14_sig = rc.d14[[2]] #4 pathways
rc.d14_sig = rc.d14_sig %>% rename(pathway = "kcat.mod")
rc.d14_sig = left_join(rc.d14_sig, ko.levels)
rc.d14_sig$leadingEdge = as.character(rc.d14_sig$leadingEdge)
#write.csv(rc.d14_sig, "rc.d14.met.gsea.csv", row.names=FALSE)

rc.d28 = kcat_gsea("DEG_vs_preABX/rc.d28.csv")
rc.d28_sig = rc.d28[[2]] #9 pathways
rc.d28_sig = rc.d28_sig %>% rename(pathway = "kcat.mod")
rc.d28_sig = left_join(rc.d28_sig, ko.levels)
rc.d28_sig$leadingEdge = as.character(rc.d28_sig$leadingEdge)
#write.csv(rc.d28_sig, "rc.d28.met.gsea.csv", row.names=FALSE)


#### now add in the day-to-day comparisons
rc.d2v4 = kcat_gsea("DEGs_byDay/rc.d2_4.csv")
rc.d2v4_sig = rc.d2v4[[2]] #15 pathways
rc.d2v4_sig = rc.d2v4_sig %>% rename(pathway = "kcat.mod")
rc.d2v4_sig = left_join(rc.d2v4_sig, ko.levels)
rc.d2v4_sig$leadingEdge = as.character(rc.d2v4_sig$leadingEdge)
#write.csv(rc.d2v4_sig, "rc.d2v4.met.gsea.csv", row.names=FALSE)

rc.d4v14 = kcat_gsea("DEGs_byDay/rc.d4_14.csv")
rc.d4v14_sig = rc.d4v14[[2]] #9 pathways
rc.d4v14_sig = rc.d4v14_sig %>% rename(pathway = "kcat.mod")
rc.d4v14_sig = left_join(rc.d4v14_sig, ko.levels)
rc.d4v14_sig$leadingEdge = as.character(rc.d4v14_sig$leadingEdge)
#write.csv(rc.d4v14_sig, "rc.d4v14.met.gsea.csv", row.names=FALSE)


rcd2=read.csv("rc.d2.met.gsea.csv")

#before this you have to open it in excel and sort it by NES
#re-level the module names to get them in the right order
rcd2$kcat.mod=factor(rcd2$kcat.mo, levels=c(rev(rcd2$kcat.mod)))

rcd2$enr = "EnrD-3"
rcd2$enr[rcd2$NES > 0] = "EnrD2"
rcd2$enr = factor(rcd2$enr, levels=c("EnrD-3","EnrD2"))
ggplot(rcd2, aes(x=kcat.mod, y=NES, fill=enr)) +
  geom_bar(stat="identity", position="dodge",width=.6) +
  geom_hline(yintercept=0,linewidth=1) +
  scale_fill_manual(values=c("darkblue","#3477eb")) +
  #scale_y_continuous(limits=c(-2.5,2.5)) +
  ylab("Normalized Enrichment Score") +
  coord_flip() +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_text(size = 12, face='bold',color='black'),
        legend.title = element_text(face='bold'))



rcd4=read.csv("rc.d4.met.gsea.csv")
rcd4$kcat.mod=factor(rcd4$kcat.mo, levels=c(rev(rcd4$kcat.mod)))

rcd4$enr = "EnrD-3"
rcd4$enr[rcd4$NES > 0] = "EnrD4"
rcd4$enr = factor(rcd4$enr, levels=c("EnrD-3","EnrD4"))
ggplot(rcd4, aes(x=kcat.mod, y=NES, fill=enr)) +
  geom_bar(stat="identity", position="dodge",width=.6) +
  geom_hline(yintercept=0,linewidth=1) +
  scale_fill_manual(values=c("darkblue","#7db7f5")) +
  #scale_y_continuous(limits=c(-2.5,2.5)) +
  ylab("Normalized Enrichment Score") +
  coord_flip() +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_text(size = 12, face='bold',color='black'),
        legend.title = element_text(face='bold'))


rcd24=read.csv("rc.d2v4.met.gsea.csv")
rcd24$kcat.mod=factor(rcd24$kcat.mo, levels=c(rev(rcd24$kcat.mod)))

rcd24$enr = "EnrD2"
rcd24$enr[rcd24$NES > 0] = "EnrD4"
rcd24$enr = factor(rcd24$enr, levels=c("EnrD2","EnrD4"))
ggplot(rcd24, aes(x=kcat.mod, y=NES, fill=enr)) +
  geom_bar(stat="identity", position="dodge",width=.6) +
  geom_hline(yintercept=0,linewidth=1) +
  scale_fill_manual(values=c("#3477eb","#7db7f5")) +
  #scale_y_continuous(limits=c(-2.5,2.5)) +
  ylab("Normalized Enrichment Score") +
  coord_flip() +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_text(size = 12, face='bold',color='black'),
        legend.title = element_text(face='bold'))


rcd414=read.csv("rc.d4v14.met.gsea.csv")
rcd414$kcat.mod=factor(rcd414$kcat.mo, levels=c(rev(rcd414$kcat.mod)))

rcd414$enr = "EnrD4"
rcd414$enr[rcd414$NES > 0] = "EnrD14"
rcd414$enr = factor(rcd414$enr, levels=c("EnrD4","EnrD14"))
ggplot(rcd414, aes(x=kcat.mod, y=NES, fill=enr)) +
  geom_bar(stat="identity", position="dodge",width=.6) +
  geom_hline(yintercept=0,linewidth=1) +
  scale_fill_manual(values=c("#a9ecf5")) +
  #scale_y_continuous(limits=c(-2.5,2.5)) +
  ylab("Normalized Enrichment Score") +
  coord_flip() +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_text(size = 12, face='bold',color='black'),
        legend.title = element_text(face='bold'))



wd.d13 = kcat_gsea("DEG_vs_preABX/wd.d13.csv")
wd.d13_sig = wd.d13[[2]] #20 pathways
wd.d13_sig = wd.d13_sig %>% rename(pathway = "kcat.mod")
wd.d13_sig = left_join(wd.d13_sig, ko.levels)
wd.d13_sig$leadingEdge = as.character(wd.d13_sig$leadingEdge)
#write.csv(wd.d13_sig, "wd.d13.met.gsea.csv", row.names=FALSE)

wd.d2 = kcat_gsea("DEG_vs_preABX/wd.d2.csv")
wd.d2_sig = wd.d2[[2]] #11 pathways
wd.d2_sig = wd.d2_sig %>% rename(pathway = "kcat.mod")
wd.d2_sig = left_join(wd.d2_sig, ko.levels)
wd.d2_sig$leadingEdge = as.character(wd.d2_sig$leadingEdge)
#write.csv(wd.d2_sig, "wd.d2.met.gsea.csv", row.names=FALSE)

wd.d14 = kcat_gsea("DEG_vs_preABX/wd.d14.csv")
wd.d14_sig = wd.d14[[2]] #10 pathways
wd.d14_sig = wd.d14_sig %>% rename(pathway = "kcat.mod")
wd.d14_sig = left_join(wd.d14_sig, ko.levels)
wd.d14_sig$leadingEdge = as.character(wd.d14_sig$leadingEdge)
#write.csv(wd.d14_sig, "wd.d14.met.gsea.csv", row.names=FALSE)

wd.d28 = kcat_gsea("DEG_vs_preABX/wd.d28.csv")
wd.d28_sig = wd.d28[[2]] #6 pathways
wd.d28_sig = wd.d28_sig %>% rename(pathway = "kcat.mod")
wd.d28_sig = left_join(wd.d28_sig, ko.levels)
wd.d28_sig$leadingEdge = as.character(wd.d28_sig$leadingEdge)
#write.csv(wd.d28_sig, "wd.d28.met.gsea.csv", row.names=FALSE)

###day by day comparisons
wd.d2_14 = kcat_gsea("DEGs_byDay/wd.d2_14.csv")
wd.d2_14_sig = wd.d2_14[[2]] #8 pathways
wd.d2_14_sig = wd.d2_14_sig %>% rename(pathway = "kcat.mod")
wd.d2_14_sig = left_join(wd.d2_14_sig, ko.levels)
wd.d2_14_sig$leadingEdge = as.character(wd.d2_14_sig$leadingEdge)
#write.csv(wd.d2_14_sig, "wd.d2_14.met.gsea.csv", row.names=FALSE)


wd.d14_28 = kcat_gsea("DEGs_byDay/wd.d14_28.csv")
wd.d14_28_sig = wd.d14_28[[2]] #8 pathways
wd.d14_28_sig = wd.d14_28_sig %>% rename(pathway = "kcat.mod")
wd.d14_28_sig = left_join(wd.d14_28_sig, ko.levels)
wd.d14_28_sig$leadingEdge = as.character(wd.d14_28_sig$leadingEdge)
#write.csv(wd.d14_28_sig, "wd.d14_28.met.gsea.csv", row.names=FALSE)


wd2=read.csv("wd.d2.met.gsea.csv")
wd2$kcat.mod=factor(wd2$kcat.mo, levels=c(rev(wd2$kcat.mod)))

wd2$enr = "EnrD-3"
wd2$enr[wd2$NES > 0] = "EnrD2"
wd2$enr = factor(wd2$enr, levels=c("EnrD-3","EnrD2"))
ggplot(wd2, aes(x=kcat.mod, y=NES, fill=enr)) +
  geom_bar(stat="identity", position="dodge",width=.6) +
  geom_hline(yintercept=0,linewidth=1) +
  scale_fill_manual(values=c("darkred","#db4237")) +
  #scale_y_continuous(limits=c(-2.5,2.5)) +
  ylab("Normalized Enrichment Score") +
  coord_flip() +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_text(size = 12, face='bold',color='black'),
        legend.title = element_text(face='bold'))



wd14=read.csv("wd.d14.met.gsea.csv")
wd14$kcat.mod=factor(wd14$kcat.mo, levels=c(rev(wd14$kcat.mod)))

wd14$enr = "EnrD-3"
wd14$enr[wd14$NES > 0] = "EnrD14"
wd14$enr = factor(wd14$enr, levels=c("EnrD-3","EnrD14"))
ggplot(wd14, aes(x=kcat.mod, y=NES, fill=enr)) +
  geom_bar(stat="identity", position="dodge",width=.6) +
  geom_hline(yintercept=0,linewidth=1) +
  scale_fill_manual(values=c("darkred","#f5958e")) +
  #scale_y_continuous(limits=c(-2.5,2.5)) +
  ylab("Normalized Enrichment Score") +
  coord_flip() +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_text(size = 12, face='bold',color='black'),
        legend.title = element_text(face='bold'))


wd28=read.csv("wd.d28.met.gsea.csv")
wd28$kcat.mod=factor(wd28$kcat.mo, levels=c(rev(wd28$kcat.mod)))

wd28$enr = "EnrD-3"
wd28$enr[wd28$NES > 0] = "EnrD28"
wd28$enr = factor(wd28$enr, levels=c("EnrD-3","EnrD28"))
ggplot(wd28, aes(x=kcat.mod, y=NES, fill=enr)) +
  geom_bar(stat="identity", position="dodge",width=.6) +
  geom_hline(yintercept=0,linewidth=1) +
  scale_fill_manual(values=c("darkred","#f5ccc9")) +
  #scale_y_continuous(limits=c(-2.5,2.5)) +
  ylab("Normalized Enrichment Score") +
  coord_flip() +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_text(size = 12, face='bold',color='black'),
        legend.title = element_text(face='bold'))




wd214=read.csv("wd.d2_14.met.gsea.csv")
wd214$kcat.mod=factor(wd214$kcat.mo, levels=c(rev(wd214$kcat.mod)))

wd214$enr = "EnrD2"
wd214$enr[wd214$NES > 0] = "EnrD14"
wd214$enr = factor(wd214$enr, levels=c("EnrD2","EnrD14"))
ggplot(wd214, aes(x=kcat.mod, y=NES, fill=enr)) +
  geom_bar(stat="identity", position="dodge",width=.6) +
  geom_hline(yintercept=0,linewidth=1) +
  scale_fill_manual(values=c("#db4237","#f5958e")) +
  #scale_y_continuous(limits=c(-2.5,2.5)) +
  ylab("Normalized Enrichment Score") +
  coord_flip() +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_text(size = 12, face='bold',color='black'),
        legend.title = element_text(face='bold'))


wd1428=read.csv("wd.d14_28.met.gsea.csv")
wd1428$kcat.mod=factor(wd1428$kcat.mo, levels=c(rev(wd1428$kcat.mod)))

wd1428$enr = "EnrD14"
wd1428$enr[wd1428$NES > 0] = "EnrD28"
wd1428$enr = factor(wd1428$enr, levels=c("EnrD14","EnrD28"))
ggplot(wd1428, aes(x=kcat.mod, y=NES, fill=enr)) +
  geom_bar(stat="identity", position="dodge",width=.6) +
  geom_hline(yintercept=0,linewidth=1) +
  scale_fill_manual(values=c("#f5958e","#f5ccc9")) +
  #scale_y_continuous(limits=c(-2.5,2.5)) +
  ylab("Normalized Enrichment Score") +
  coord_flip() +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_text(size = 12, face='bold',color='black'),
        legend.title = element_text(face='bold'))







