library(tidyverse)

setwd("/DIRECTORY/")

data = read.csv("a1_a3_allCFUs.csv")

#change 0s to 1000 for plotting
data$meanCFUpg_BHIS[which(data$meanCFUpg_BHIS == 0)] = 1e+3

#subset data for first 14 days after abx
d14 = subset(data, Rec_day_adj <= 14)

#subset the female cohorts (A1C1, A3C123)
d14_c123 = subset(d14, Cohort != 5 & Cohort != 6) %>% subset(Cohort != 4)
d14_c123 = subset(d14_c123, Rec_day_adj != 13) #only one cohort has a D13, so just get rid of that day

d14_c123 %>% group_by(Trt_noInf, Rec_day_adj) %>% tally()
#for females, n=10-24/group

#subset the male cohorts (A3C4)
c4 = subset(d14, Cohort == 4)
c4 %>% group_by(Trt_noInf, Rec_day_adj) %>% tally()
#for males, n=6/group


#calculate mean CFUs/g for each treatment/timepoint for female cohort
d14.summ = d14_c123 %>% group_by(Trt_noInf, Rec_day_adj, Diet, Abx) %>% 
  summarise(
    mean = mean(meanCFUpg_BHIS,na.rm=TRUE),
    sd = sd(meanCFUpg_BHIS, na.rm=TRUE))

p = ggplot(d14.summ) + aes(x = Rec_day_adj, y = mean, col = Trt_noInf) +
  geom_errorbar(aes(ymax=mean+sd,ymin=mean-sd),
                position=position_dodge(.5),
                width=2, linewidth=.8) +
  scale_y_log10() +
  scale_y_log10(breaks=c(1e+03,1e+04,1e+05,1e+06,1e+07,1e+08,1e+09,1e+10,1e+11, 1e+12),limits=c(1e+03,1e+12)) + 
  geom_line(alpha=0.8) + geom_jitter(position=position_dodge(0.5)) +
  ylab("CFU/g fecal material, BHIS") + theme_bw() +
  scale_x_continuous(breaks=c(-3,0,4,7,10,14),limits=c(-4,15)) +
  scale_color_manual(values=c("#00A9FF",  "#7CAE00",
                              "#F8766D", "#CD9600")) +
  facet_grid(Diet~Abx) +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  theme(legend.title = element_text(size = 9, face="bold"), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"),
        axis.title.x = element_blank(),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p

#ggsave(file="FILE.svg", plot=p, width=5, height=4)

#calculate mean CFU/g for each treatment/timepoint for male cohort
c4.summ = c4 %>% group_by(Treatment, Rec_day_adj, Diet, Abx) %>% 
  summarise(
    mean = mean(meanCFUpg_BHIS),
    sd = sd(meanCFUpg_BHIS, na.rm=TRUE))

p1 = ggplot(c4.summ) + aes(x = Rec_day_adj, y = mean, col = Treatment) +
  geom_errorbar(aes(ymax=mean+sd,ymin=mean-sd),
                position=position_dodge(.5),
                width=2, linewidth=.8) +
  scale_y_log10(breaks=c(1e+03,1e+04,1e+05,1e+06,1e+07,1e+08,1e+09,1e+10,1e+11, 1e+12),limits=c(1e+03,1e+12)) + geom_line() + geom_point() +
  ylab("CFU/g fecal material, BHIS") + theme_bw() +
  scale_x_continuous(breaks=c(-3,0,4,7,10,14),limits=c(-4,15)) +
  scale_color_manual(values=c("#00A9FF","#F8766D")) +
  facet_grid(Diet~Abx) +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  theme(legend.title = element_text(size = 9, face="bold"), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"),
        axis.title.x = element_blank(),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
p1

#ggsave(file="FILE.svg", plot=p1, width=3.5, height=4)