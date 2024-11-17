library(tidyverse)

setwd("~/Documents/Chang-Bergelson Lab/HF Diet Abx/Aim 1/Metagenomics/comm_analysis/")

#this script counts up all the overlaps in the DEGs identified in 05.deseq2_KO_filt.R
#which were used to manually create Figure S2D, S2E, S2G, and S2H in Inkscape

rc2 = read.csv("DEG_vs_preABX/rc.d2.DEG.csv")
rc2.dep = rc2$KO[rc2$log2FoldChange < 0] #i.e. depleted at d2 rel to pre-abx
rc2.enr = rc2$KO[rc2$log2FoldChange > 0] #i.e. enriched at d2 rel to pre-abx

#depleted/enriched at D4 rel to pre-abx on RC
rc4 = read.csv("DEG_vs_preABX/rc.d4.DEG.csv")
rc4.dep = rc4$KO[rc4$log2FoldChange < 0]
rc4.enr = rc4$KO[rc4$log2FoldChange > 0]

#depleted at D14 rel to pre-abx
rc14 = read.csv("DEG_vs_preABX/rc.d14.DEG.csv")
rc14.dep = rc14$KO[rc14$log2FoldChange < 0]

#depleted at D28 rel to pre-abx
rc28 = read.csv("DEG_vs_preABX/rc.d28.DEG.csv")
rc28.dep = rc28$KO[rc28$log2FoldChange < 0]

#depleted at D2, still depleted at D4
rc.2_4.dep = rc4.dep[rc4.dep %in% rc2.dep] 
#438/658 KOs depleted at D2 are still missing at D4

#depleted at D2, still depleted at D4, still depleted at D14
rc.2_4_14.dep = rc14.dep[rc14.dep %in% rc.2_4.dep] 
#all 6 d14 KOs are a subset of the ones that were gone at d2 and d4

#depleted from D2 through D28
rc.2_4_14_28.dep = rc28.dep[rc28.dep %in% rc.2_4_14.dep]
#only 1 KO

rc28.dep[rc28.dep %in% rc14.dep]
rc28.dep[rc28.dep %in% rc4.dep] #the other KO depleted at D28 was also depleted at D4
rc28.dep[rc28.dep %in% rc2.dep]




#depleted/enriched at d2 rel to pre-abx on WD
wd2 = read.csv("DEG_vs_preABX/wd.d2.DEG.csv")
wd2.dep = wd2$KO[wd2$log2FoldChange < 0]
wd2.enr = wd2$KO[wd2$log2FoldChange > 0]

#depleted/enriched at d14 rel to pre-abx
wd14 = read.csv("DEG_vs_preABX/wd.d14.DEG.csv")
wd14.dep = wd14$KO[wd14$log2FoldChange < 0]
wd14.enr = wd14$KO[wd14$log2FoldChange > 0]

#depleted/enriched at d28 rel to pre-abx
wd28 = read.csv("DEG_vs_preABX/wd.d28.DEG.csv")
wd28.dep = wd28$KO[wd28$log2FoldChange < 0]
wd28.enr = wd28$KO[wd28$log2FoldChange > 0]

#nKOs depleted at d2 that were still depleted at d14
wd.2_14.dep = wd14.dep[wd14.dep %in% wd2.dep] 
#641/939 genes depleted at D2 are still missing at D14
#3 new KOs are depleted at D14 that weren't depleted at d2

#KOs depleted at D2 and no other timepoints
wd2.dep.unique = wd2.dep[!wd2.dep %in% wd14.dep]
wd2.dep.unique = wd2.dep.unique[!wd2.dep.unique %in% wd28.dep]
#write.csv(wd2.dep.unique, "wd.rec.2_14.csv", row.names=FALSE)

#KOs depleted from D2 through D28
wd.2_14_28.dep = wd28.dep[wd28.dep %in% wd.2_14.dep] 
#273/641 KOs depleted at both d2 and d14 are still missing at d28
#write.csv(wd28.dep, "wd.d28.dep.csv", row.names=FALSE)

#KOs depleted at D14 and no other timepoints
wd14.dep.unique = wd14.dep[!wd14.dep %in% wd2.dep] 
wd14.dep.unique = wd14.dep.unique[!wd14.dep.unique %in% wd28.dep]
#only 3 KOs

#KOs depleted at D14 that stayed depleted to D28
wd.14_28.dep = wd28.dep[wd28.dep %in% wd14.dep] 
#all of these were also depleted at d2

wd.2_28.dep = wd28.dep[wd28.dep %in% wd2.dep] 
#all of the D28 KOs were depleted at D2


wd.2_14.enr = wd14.enr[wd14.enr %in% wd2.enr] 
#98 D14 KOs enriched at D2 and D14

wd2.unique = wd2.enr[!wd2.enr %in% wd14.enr] 
wd2.unique = wd2.unique[!wd2.unique %in% wd28.enr]
#69 KOs enriched at only D2 on WD

wd.2_14_28.enr = wd28.enr[wd28.enr %in% wd.2_14.enr] 
#10 D28 KOs enriched from D2 throuhg D28

wd.14_28.enr = wd28.enr[wd28.enr %in% wd14.enr] 
#36 KOs enriched at D14
#26 unique to D14, 10 common to D2 and D14

wd14.unique = wd14.enr[!wd14.enr %in% wd2.enr] 
wd14.unique = wd14.unique[!wd14.unique %in% wd28.enr]
#26 KOs uniquely enriched at D14

wd.2_28.enr = wd28.enr[wd28.enr %in% wd2.enr] 
#14 enriched at D2
#4 unique to 2, 10 common to 2 and 14

#to figure out if there are any unique to D28... 
wd28.unique = wd28.enr[!wd28.enr %in% wd2.enr] 
wd28.unique = wd28.unique[!wd28.unique %in% wd14.enr]
#6 unique to D28


wd.rec.14_28 = wd2.dep[wd2.dep %in% wd14.dep]
wd.rec.14_28 = wd14.dep[!wd14.dep %in% wd28.dep]
#write.csv(wd.rec.14_28, "wd.rec.14_28.csv",row.names=FALSE)




#KOs that are depleted at D2 but recover by D4
rc.d2.unq = rc2.dep[!rc2.dep %in% rc4.dep]
rc.d2.unq = rc.d2.unq[!rc.d2.unq %in% rc14.dep]
rc.d2.unq = rc.d2.unq[!rc.d2.unq %in% rc28.dep]
rc.rec.2_4 = rc.d2.unq 
#write.csv(rc.rec.2_4, "rc.rec.2_4.csv", row.names=FALSE)

#KOs depleted from D2 through D4 but recover by D14
rc.rec.4_14 = rc.2_4.dep[!rc.2_4.dep %in% rc14.dep]
rc.rec.4_14 = rc.rec.4_14[!rc.rec.4_14 %in% rc28.dep] #rec'd bw D4 and D14
#write.csv(rc.rec.4_14, "rc.rec.4_14.csv", row.names=FALSE)

#KOs depleted through D14 but recover by D28
rc.rec.14_28 = rc.2_4_14.dep[!rc.2_4_14.dep %in% rc28.dep] #5 KO
#write.csv(rc.rec.14_28, "rc.rec.14_28.csv", row.names=FALSE)


#the unique D4 depleted KOs that recover by D14 are...
rc.d4.unq = rc4.dep[!rc4.dep %in% rc2.dep] #these are the ones newly depleted at D4
#they also all recover by D14
rc.rec.all4 = rc.d4.unq[!rc.d4.unq %in% rc28.dep]
#write.csv(rc.rec.all4, "rc.rec.4unq.csv", row.names=FALSE)














