library(tidyverse)

setwd("~/Documents/Chang-Bergelson Lab/HF Diet Abx/Aim 3 col resistance/Data/")

water = read.csv("H2O.csv")
water = dplyr::select(water, -Date)
water.wide = spread(water, key="Desc",value="Water.weight")


water.wide$deltWater = water.wide$pre-water.wide$post

water.wide$water.p.mouse.p.day = water.wide$deltWater/water.wide$nMice/3

water.wide$trt = paste0(water.wide$Diet,"-", water.wide$Abx)

#Figure S1A
p = ggplot(water.wide, aes(x=Diet,y=water.p.mouse.p.day, fill=trt,color=trt)) + 
  geom_boxplot() + #facet_grid(.~Cohort) + 
  ylab("Water Consumption/Mouse/Day (ml)") +
  scale_fill_manual(values=c("blue","white","red","white")) +
  scale_color_manual(values=c("black","blue","black","red")) +
  theme(panel.grid.major = element_line(color="gray90"),
        #panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_text(size = 9, face="bold"), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.title.y = element_blank(),
        panel.spacing.x = unit(2, "lines"),
        strip.text = element_blank(),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p

#ggsave(file="water.svg", plot=p, width=3, height=2.5)

#stats done in prism

#write.csv(water.wide, "~/Documents/WD_ABX manuscript/nature_drafts/July_alt/supp figs/supp_16s/water_data.csv",row.names=FALSE)
