library(tidyverse)
library(RColorBrewer)

#this script calculates and plots the NUMBER of sig. depleted KOs mapping to 
#each KEGG family across timepoints and diets

setwd("/DIRECTORY/")

ko = read.csv("ko.mod4.csv")
kfam = dplyr::select(ko, c("KO","kfam.mod"))
kfam = unique(kfam)
kcat = dplyr::select(ko, c("KO","kcat.mod"))
kcat = unique(kcat)

#import lists of KOs
rc24rec = read.csv("DEGs_byDay/rc.rec.2_4.csv") #depleted at d2 but recover by d4
rc414rec = read.csv("DEGs_byDay/rc.rec.4_14.csv") #depleted from d2-d4, but recovery by d14
rc1428rec = read.csv("DEGs_byDay/rc.rec.14_28.csv") #depleted throuhg d14, recover by d28
rc4rec = read.csv("DEGs_byDay/rc.rec.4unq.csv") #newly depleted at d4, recover by d14

#map to KEGG family
rc24rec.kfam = left_join(rc24rec, kfam, multiple="all")
rc414rec.kfam = left_join(rc414rec, kfam, multiple="all")
rc1428rec.kfam = left_join(rc1428rec, kfam, multiple="all")
rc4rec.kfam = left_join(rc4rec, kfam, multiple="all")

#filter out ones without informative KEGG mapping
#tally how many KOs per KEGG family
rc24.kfam.tally = rc24rec.kfam %>% group_by(kfam.mod) %>% 
  subset(kfam.mod != "NA") %>%
  subset(kfam.mod != "09194 Poorly characterized") %>%
  tally()
rc24.kfam.tally$day = "D2"

rc414.kfam.tally = rc414rec.kfam %>% group_by(kfam.mod) %>% 
  subset(kfam.mod != "NA") %>%
  subset(kfam.mod != "09194 Poorly characterized") %>%
  tally()
rc414.kfam.tally$day = "D2/4"

rc4.kfam.tally = rc4rec.kfam %>% group_by(kfam.mod) %>% 
  subset(kfam.mod != "NA") %>%
  subset(kfam.mod != "09194 Poorly characterized") %>%
  tally()
rc4.kfam.tally$day = "D4"

rc1428.kfam.tally = rc1428rec.kfam %>% group_by(kfam.mod) %>% 
  subset(kfam.mod != "NA") %>%
  subset(kfam.mod != "09194 Poorly characterized") %>%
  tally()
rc1428.kfam.tally$day = "D14"

#combine the tallies
rc.all_rec = rbind(rc24.kfam.tally, rc414.kfam.tally, rc4.kfam.tally)
rc.all_rec$day = factor(rc.all_rec$day, levels=c("D2","D2/4","D4"))

#map to KEGG system for ordering reasons
ksys = dplyr::select(ko, c("kfam.mod", "ksys.mod"))
ksys = unique(ksys)
ksys$ksys.mod = factor(ksys$ksys.mod, levels=c("Metabolism","Cellular processes",
                                               "Genetic information processing",
                                               "Environmental information processing",
                                               "Not Included in Pathway or Brite"))
rc.all_rec.ko = left_join(rc.all_rec, ksys)

#order them by KEGG system
#write.csv(rc.all_rec.ko, "rc.all_rec.ko.csv",row.names=FALSE)
#rc.all_rec.ko2 = read.csv("rc.all_rec.ko.csv")
rc.all_rec.ko2 = rc.all_rec.ko

rc.all_rec.ko2$kfam.mod <- factor(rc.all_rec.ko2$kfam.mod,
                               levels = unique(rc.all_rec.ko2$kfam.mod), ordered = TRUE)

cols=c("#7FC97F","#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#ff9900", "grey70",  "#1B9E77",
       "#D95F02", "#f78b83","#7570B3",  "#66A61E", "#E6AB02", "#A6761D","#666666", "#A6CEE3", "#1F78B4",
       "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "turquoise3", "#FF7F00", "#7FC97F", "#BEAED4", "#FDC086",
       "#FFFF99")

#plot the subsets of KOs depleted at each timepoint colored by KEGG family
#This is Figure S2I
p1 = ggplot(data = rc.all_rec.ko2, aes(x=day, fill=kfam.mod, weight=n)) + 
  geom_bar(position="stack") + 
  scale_fill_manual(values=cols) +
  scale_y_continuous(breaks=c(seq(0,500,100)),limits=c(0,530)) +
  theme(#axis.line = element_line(colour = "black"),
    panel.grid.major = element_line(color="lightgrey"),
    panel.grid.minor = element_line(color="lightgrey"),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        legend.title = element_text(face='bold'),
        legend.position="none") +
  guides(fill=guide_legend(title="KEGG Family"),ncol=1)
p1
#ggsave(file="RC_depleted_kfam.svg", plot=p1, width=2.5, height=3)



wd214rec = read.csv("DEGs_byDay/wd.rec.2_14.csv") #depleted at d2 but recovered by d14
wd1428rec = read.csv("DEGs_byDay/wd.rec.14_28.csv") #depleted through d14, recovered by d28
wd28dep = read.csv("DEGs_byDay/wd.d28.dep.csv") #depleted through d28

#map to KEGG family
wd214rec.kfam = left_join(wd214rec, kfam, multiple="all")
wd1428rec.kfam = left_join(wd1428rec, kfam, multiple="all")
wd28dep.kfam = left_join(wd28dep, kfam, multiple="all")

#fiter out uninformative KEGG mappings
#tally nKOs per KEGG family
wd214.kfam.tally = wd214rec.kfam %>% group_by(kfam.mod) %>% 
  subset(kfam.mod != "NA") %>%
  subset(kfam.mod != "09194 Poorly characterized") %>%
  tally()
wd214.kfam.tally$day = "D2/4-D14"

wd1428.kfam.tally = wd1428rec.kfam %>% group_by(kfam.mod) %>% 
  subset(kfam.mod != "NA") %>%
  subset(kfam.mod != "09194 Poorly characterized") %>%
  tally()
wd1428.kfam.tally$day = "D14-D28"

wd28dep.kfam.tally = wd28dep.kfam %>% group_by(kfam.mod) %>% 
  subset(kfam.mod != "NA") %>%
  subset(kfam.mod != "09194 Poorly characterized") %>%
  tally()
wd28dep.kfam.tally$day = "D28-dep"

#combine tallies
wd.all_rec = rbind(wd214.kfam.tally, wd1428.kfam.tally, wd28dep.kfam.tally)
wd.all_rec$day = factor(wd.all_rec$day, levels=c("D2/4-D14","D14-D28", "D28-dep"))

#map to KEGG system for ordering
ksys = dplyr::select(ko, c("kfam.mod", "ksys.mod"))
ksys = unique(ksys)
ksys$ksys.mod = factor(ksys$ksys.mod, levels=c("Metabolism","Cellular processes",
                                               "Genetic information processing",
                                               "Environmental information processing",
                                               "Viral protein families",
                                               "Not Included in Pathway or Brite"))
wd.all_rec.ko = left_join(wd.all_rec, ksys)

#order them by KEGG system from excel so they match
#write.csv(wd.all_rec.ko, "wd.all_rec.ko.csv",row.names=FALSE)
wd.all_rec.ko2 = read.csv("wd.all_rec.ko.csv")

wd.all_rec.ko2$kfam.mod <- factor(wd.all_rec.ko2$kfam.mod,
                               levels = unique(wd.all_rec.ko2$kfam.mod), ordered = TRUE)

wd.all_rec.ko2$day = factor(wd.all_rec.ko2$day, levels=c("D2/4-D14","D14-D28", "D28-dep"))

wd.cols = c(cols[1:22],"#FFFF99",cols[23:27])

#plotting WD KOs depleted at each timepoint
#This is Figure S2L
p2 = ggplot(data = wd.all_rec.ko2, aes(x=day, fill=kfam.mod, weight=n)) + 
  geom_bar(position="stack") + 
  scale_fill_manual(values=wd.cols) +
  scale_y_continuous(breaks=c(seq(0,500,100)),limits=c(0,530)) +
  theme(panel.grid.major = element_line(color="lightgrey"),
    panel.grid.minor = element_line(color="lightgrey"),
    panel.border = element_rect(linewidth = 1,fill=NA),
    panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        legend.title = element_text(face='bold'),
        legend.text  = element_text(size = 7),
        legend.key.size = unit(0.6, "lines")) +
  guides(fill=guide_legend(title="KEGG Family",ncol=1))
p2
#ggsave(file="WD_depleted_kfam.svg", plot=p2, width=4.9, height=3)





