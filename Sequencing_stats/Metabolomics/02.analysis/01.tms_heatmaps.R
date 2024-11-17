library(tidyverse)
library(scales)

setwd("DIRECTORY/")

#TMS panel
gcms = read.csv("TMS.csv", check.names = FALSE)

#set the NAs to 0
gcms[is.na(gcms)] = 0

#convert table to tidy/long form
gcms.tidy = gather(gcms, key="Compound",value="value", -SampleName)

#this file categorizes each compound in broader groups
class = read.csv("compound_classes.csv")
gcms.tidy = left_join(gcms.tidy, class)

#import metadata
meta = read.csv("metabolomics_meta.csv")
meta = rename(meta, "SampleName" = Sample.ID)

#merge abundance table with metadata
gcms.meta = left_join(gcms.tidy, meta)

gcms.meta$Class = as.factor(gcms.meta$Class)
gcms.meta$Rec_day_adj = factor(gcms.meta$Rec_day_adj, levels=c("-3","0","3","5","7","11","14","28"))

#make a dummy variable (mouse-compound) for comparing each mouse to its baseline
gcms.meta1 = gcms.meta %>% mutate("normID" = paste(Mouse, Compound, sep="_"))

#there are cecal samples in this dataset that we are not interested in
#so filter to just fecal
fecal = subset(gcms.meta1, Type=="fecal")
fecal.meta = subset(meta, Type=="fecal")

#each mouse-compound has a unique pre-abx value. store these, then normalize

#1946 and 1953 are the only fecals without baselines
#can either normalize to the average baseline for each treatment group
#or can just do the averaged metrics

#fecal, individual mouse results

#filter out Mouse 1946 and 1953 because they did not have baseline pre-ABX values
fecal.complete = subset(fecal, Mouse!=1946 & Mouse!=1953)

#gather the pre-ABX values
fecal.pre = subset(fecal.complete, Rec_day_adj=="-3")
#pick out just the mouse-compound ID and the pre-ABX values
fecal.pre = data.frame("baseline" = fecal.pre$value, "normID" = fecal.pre$normID)

#merge the baseline values for each compound/mouse to the other values
fecal.bl = left_join(fecal.complete, fecal.pre)

###Normalize to pre-abx baseline and calculate log2foldchange 
#set the 0s to a very small number for log normalizing reasons
fecal.bl1 = fecal.bl %>% mutate("log2foldchange" = log2((value + 0.01) / (baseline  +0.01)))

#order the compounds by class
fecal.bl1 = fecal.bl1 %>% 
  arrange(desc(Class)) %>% 
  mutate(Compound = fct_inorder(factor(Compound, ordered=TRUE)))

#write.csv(fecal.bl1, "norm_fecal_TMS.csv",row.names=FALSE)

#subset by diet for visualization
rc = subset(fecal.bl1, Treatment=="RC-ABX" & Rec_day_adj != -3)
wd = subset(fecal.bl1, Treatment=="WD-ABX" & Rec_day_adj != -3)

#re-order the compounds so they're grouped by class
rc = rc %>% 
  arrange(desc(Class), desc(Compound)) %>% 
  mutate(Compound = fct_inorder(factor(Compound, ordered=TRUE)))

wd = wd %>% 
  arrange(desc(Class), desc(Compound)) %>% 
  mutate(Compound = fct_inorder(factor(Compound, ordered=TRUE)))


####restrict to the compounds present in the first metabolomics run
# (the second run added new compounds)
rc.run1 = subset(rc, Run1=="y")
#re-order the compounds so they're grouped by class
rc.run1 = rc.run1 %>% 
  arrange(desc(Class), desc(Compound)) %>% 
  mutate(Compound = fct_inorder(factor(Compound, ordered=TRUE)))

#Figure S3A
rc.plot = ggplot(rc.run1, aes(x=as.factor(Mouse), y = Compound, fill=log2foldchange)) + 
  geom_tile(color="grey") + facet_wrap(Rec_day_adj~., ncol=6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(legend.title = element_text(size = 7), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.4, "lines")) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        panel.spacing.x = unit(1, "lines"),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=8,angle=45,vjust=1,hjust=1),
        legend.title = element_text(face='bold'),
        strip.text.x = element_text(face="bold",size=12),
        strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
rc.plot
#ggsave(file="RC_heatmap_byMouse.svg", plot=rc.plot, width=10, height=8)


wd.run1 = subset(wd, Run1=="y")
wd.run1 = wd.run1 %>% 
  arrange(desc(Class), desc(Compound)) %>% 
  mutate(Compound = fct_inorder(factor(Compound, ordered=TRUE)))

#Figure S3A
wd.plot = ggplot(wd.run1, aes(x=as.factor(Mouse), y = Compound, fill=log2foldchange)) + 
  geom_tile(color="grey") + facet_wrap(Rec_day_adj~., ncol=6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(legend.title = element_text(size = 7), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.4, "lines")) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        panel.spacing.x = unit(1, "lines"),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=8,angle=45,vjust=1,hjust=1),
        legend.title = element_text(face='bold'),
        strip.text.x = element_text(face="bold",size=12),
        strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
wd.plot
#ggsave(file="WD_heatmap_byMouse.svg", plot=wd.plot, width=6, height=8)




#now average these values across mice

#restrict to compounds from first metabolomics run
gcms.run1 = subset(fecal.bl1, Run1=="y")

#average across mice within timepoint/diet/compound
gcms.summ = gcms.run1 %>% group_by(Rec_day_adj, Treatment, Compound) %>% 
  summarise(
    mean_value = mean(value),
    sd = sd(value, na.rm=TRUE))

#okay, and now calculate fold change relative to mean baseline and log2 it
#dummy variable for compound-treatment for matching up to baseline
gcms.summ = mutate(gcms.summ, normID = paste(Compound, Treatment, sep="_" ))

#pre-ABX baseline values for each compound
pre = subset(gcms.summ, Rec_day_adj=="-3")
#just take the pre-ABX value and the compound-treatment ID for matching up
pre = data.frame("baseline" = pre$mean_value, "normID" = pre$normID)

#match up pre-ABX values with later timepoints
gcms.summ1 = left_join(gcms.summ, pre)

#calculate log2 fold change
gcms.summ1 = gcms.summ1 %>% mutate("log2foldchange" = log2((mean_value + 0.1) / (baseline  +0.1)))

#merge w compound class for ordering
gcms.summ1 = left_join(gcms.summ1, class)

#re-order the compounds so they're grouped by class
gcms.summ2 = gcms.summ1 %>% 
  arrange(desc(Class), desc(Compound)) %>% 
  mutate(Compound = fct_inorder(factor(Compound, ordered=TRUE)))

#filter out the pre-ABX timepoint since this normalizes to itself
fecal.mean = subset(gcms.summ2, Rec_day_adj != -3)

#subset RC
rc.fecal = subset(fecal.mean, Treatment=="RC-ABX")

#write.csv(rc.fecal, "rc.fecal.run1.csv",row.names=FALSE)
#rc.fecal = read.csv("rc.fecal.run1.csv")
rc.fecal$Rec_day_adj = factor(rc.fecal$Rec_day_adj, levels=c("0","3","5","7","11","14"))
#order by compound class
rc.fecal$Compound <- factor(rc.fecal$Compound, levels=unique(rc.fecal$Compound),ordered = TRUE)

#Figure 2C
p = ggplot(rc.fecal, aes(x=Rec_day_adj, y = Compound, fill=log2foldchange)) + 
  geom_tile(color="grey") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits=c(-11,11)) + 
  theme(legend.title = element_text(size = 7), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.4, "lines")) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        panel.spacing.x = unit(1, "lines"),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        legend.title = element_text(face='bold'),
        strip.text.x = element_text(face="bold",size=12),
        strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p
#ggsave(file="metab_heatmap_RC.svg", plot=p, width=6, height=8)

#subset WD
wd.fecal = subset(fecal.mean, Treatment=="WD-ABX")

#write.csv(wd.fecal, "wd.fecal.csv",row.names=FALSE)
#wd.fecal = read.csv("wd.fecal.csv")
wd.fecal$Rec_day_adj = factor(wd.fecal$Rec_day_adj, levels=c("0","7","11","14"))
wd.fecal$Compound <- factor(wd.fecal$Compound, levels=unique(wd.fecal$Compound),ordered = TRUE)

#Figure 2C
p1 = ggplot(wd.fecal, aes(x=Rec_day_adj, y = Compound, fill=log2foldchange)) + 
  geom_tile(color="grey") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",limits=c(-11,11)) + 
  theme(legend.title = element_text(size = 7), 
        legend.text  = element_text(size = 5),
        legend.key.size = unit(0.4, "lines")) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        panel.spacing.x = unit(1, "lines"),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=10),
        legend.title = element_text(face='bold'),
        legend.position = "none",
        strip.text.x = element_text(face="bold",size=12),
        strip.background.x = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p1
#ggsave(file="metab_heatmap_WD.svg", plot=p1, width=2.5, height=8)



