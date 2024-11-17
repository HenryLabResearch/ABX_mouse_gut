library(tidyverse)
library(DescTools)

setwd("DIRECTORY/")

lb = read.csv("a3_inf.csv")
lb$meanCFUpg_LB = as.numeric(lb$meanCFUpg_LB)


#type.adj2 is there for cohorts where we didn't take actual fecal samples on t=96,
#we only have the sac samples.
#so the idea was, for cohorts where we have fecal, use that, for cohorts where we
#only have col-L, use that as fecal proxy
fecal = subset(lb, Type.adj2=="fecal")
fecal = fecal[!is.na(fecal$meanCFUpg_LB),]

#write.csv(fecal,"St_inf_fecal.csv",row.names=FALSE)

t96 = subset(lb, hpi==96)
t96 = subset(t96, Type!="fecal")

t96$Type = factor(t96$Type, levels=c("ile-L","cec-L","col-L","liv","MLN","spl"))

#Figure S8E
p = ggplot(t96) + aes(x = Treatment, y = (meanCFUpg_LB+300), fill = Treatment) + 
  facet_wrap(Type~.) +
  geom_boxplot() + #geom_point() + 
  scale_y_log10(limits=c(3e+02,2e+11),breaks=c(1e+03,1e+04,1e+05,1e+06,1e+07,1e+08,1e+09,1e+10,1e+11)) +
  scale_fill_manual(values=c("blue","#5dbde5","red","#ff9790")) +
  ylab("CFU/g fecal material") + xlab("Hours post-infection") +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),
                     panel.background = element_blank()) +
  theme(plot.title = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1),
        legend.position="none",
        legend.title = element_text(face='bold'))
p
#ggsave(file="st_cfu_tissues.svg", plot=p, width=5, height=4)
#write.csv(t96,"t96_tissues.csv",row.names=FALSE)


#Fecal St AUC
st = read.csv("St_inf_fecal.csv")

auc = st %>% group_by(Mouse, Treatment) %>% 
  summarize(auc = AUC(x=hpi, y=meanCFUpg_LB))

#write.csv(auc, "St_auc_c468.csv",row.names=FALSE)
#analyze in prism


