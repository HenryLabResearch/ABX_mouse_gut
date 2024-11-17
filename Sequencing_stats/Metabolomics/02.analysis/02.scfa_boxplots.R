library(tidyverse)
library(scales)

setwd("DIRECTORY/")

scfa.quant = read.csv("SCFA.csv", check.names = FALSE)

scfa.class = read.csv("scfa_compound_classes.csv", check.names=FALSE)

#convert abundance data to long form
scfa.q.tidy = gather(scfa.quant, key="Compound",value="value", -sampleid)

#import metadata
meta = read.csv("metabolomics_meta.csv")
meta = rename(meta, "sampleid" = Sample.ID)

#merge tidy data to metadata
scfa.q.meta = left_join(scfa.q.tidy, meta)

#set recovery day as factor
scfa.q.meta$Rec_day_adj = factor(scfa.q.meta$Rec_day_adj, levels=c("-3","14","28"))

#analyze only cecal SCFA abundances
cecal = subset(scfa.q.meta, Type=="cecal")

#subset each compound and plot
acetate = subset(cecal, Compound=="acetate")
#Figure 2D
pA = ggplot(acetate, aes(x=Rec_day_adj, y=value, fill=Treatment)) +
  geom_boxplot() + ylab("Acetate Concentration (mM)") +
  xlab("Exp. Day") + facet_wrap(Treatment~.) +
  theme_bw() + scale_fill_manual(values=c("cornflowerblue","#d6644b")) +
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size = 9, face="bold"), 
        legend.text  = element_text(size = 8),
        #legend.key.size = unit(0.8, "lines"),
        axis.title=element_blank(),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
pA
#ggsave(file="acetate_comb.svg", plot=pA, width=4, height=2)



but = subset(cecal, Compound=="butyrate")
#Figure 2E
pB = ggplot(but, aes(x=Rec_day_adj, y=value, fill=Treatment)) +
  geom_boxplot() + ylab("Butyrate Concentration (mM)") +
  xlab("Exp. Day") + facet_wrap(Treatment~.) +
  theme_bw() + scale_fill_manual(values=c("cornflowerblue","#d6644b")) +
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size = 9, face="bold"), 
        legend.text  = element_text(size = 8),
        #legend.key.size = unit(0.8, "lines"),
        axis.title=element_blank(),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
pB
#ggsave(file="butyrate_comb.svg", plot=pB, width=4, height=2)


prop = subset(cecal, Compound=="propionate")
#Figure 2F
pP = ggplot(prop, aes(x=Rec_day_adj, y=value, fill=Treatment)) +
  geom_boxplot() + ylab("Propionate Concentration (mM)") +
  xlab("Exp. Day") + facet_wrap(Treatment~.) +
  theme_bw() + scale_fill_manual(values=c("cornflowerblue","#d6644b")) +
  theme(panel.background = element_blank()) +
  theme(legend.title = element_text(size = 9, face="bold"), 
        legend.text  = element_text(size = 8),
        #legend.key.size = unit(0.8, "lines"),
        axis.title=element_blank(),
        strip.text = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white",color="black",linewidth=1),
        panel.border = element_rect(color="black",fill="NA",linewidth=1))
pP
#ggsave(file="propionate_comb.svg", plot=pP, width=4, height=2)


#statistics on differential abundance were calculated in prism

### Statistics on differential abundance of metabolites over time

ace.rc = subset(acetate, Treatment=="RC-ABX")
ace.rc.3 = subset(ace.rc, Rec_day_adj=="-3")
ace.rc.3 = ace.rc.3$value

ace.rc.14 = subset(ace.rc, Rec_day_adj=="14")
ace.rc.14 = ace.rc.14$value

ace.rc.28 = subset(ace.rc, Rec_day_adj=="28")
ace.rc.28 = ace.rc.28$value

ace.wd = subset(acetate, Treatment=="WD-ABX")
ace.wd.3 = subset(ace.wd, Rec_day_adj=="-3")
ace.wd.3 = ace.wd.3$value

ace.wd.14 = subset(ace.wd, Rec_day_adj=="14")
ace.wd.14 = ace.wd.14$value

ace.wd.28 = subset(ace.wd, Rec_day_adj=="28")
ace.wd.28 = ace.wd.28$value

shapiro.test(ace.rc.3) #normal
shapiro.test(ace.rc.14) #normal
shapiro.test(ace.rc.28) #normal

var.test(x=ace.rc.3,y=ace.rc.14) #equal variance
var.test(x=ace.rc.3,y=ace.rc.28) #equal variance

t.test(ace.rc.3,ace.rc.14, alternative="two.sided", var.equal=TRUE) #p=.5, NS
t.test(ace.rc.3,ace.rc.28, alternative="two.sided", var.equal=TRUE) #p=0.36, NS

shapiro.test(ace.wd.3) #normal
shapiro.test(ace.wd.14) #normal
shapiro.test(ace.wd.28) #normal

var.test(x=ace.wd.3,y=ace.wd.14) #equal variance
var.test(x=ace.wd.3,y=ace.wd.28) #equal variance

t.test(ace.wd.3,ace.wd.14, alternative="two.sided", var.equal=TRUE) #p=8.89e-05
t.test(ace.wd.3,ace.wd.28, alternative="two.sided", var.equal=TRUE) #p=0.007325


### Butyrate now

but.rc = subset(but, Treatment=="RC-ABX")
but.rc.3 = subset(but.rc, Rec_day_adj=="-3")
but.rc.3 = but.rc.3$value

but.rc.14 = subset(but.rc, Rec_day_adj=="14")
but.rc.14 = but.rc.14$value

but.rc.28 = subset(but.rc, Rec_day_adj=="28")
but.rc.28 = but.rc.28$value

but.wd = subset(but, Treatment=="WD-ABX")
but.wd.3 = subset(but.wd, Rec_day_adj=="-3")
but.wd.3 = but.wd.3$value

but.wd.14 = subset(but.wd, Rec_day_adj=="14")
but.wd.14 = but.wd.14$value

but.wd.28 = subset(but.wd, Rec_day_adj=="28")
but.wd.28 = but.wd.28$value

shapiro.test(but.rc.3) #normal
shapiro.test(but.rc.14) #normal
shapiro.test(but.rc.28) #normal

var.test(x=but.rc.3,y=but.rc.14) #equal variance
var.test(x=but.rc.3,y=but.rc.28) #equal variance

t.test(but.rc.3,but.rc.14, alternative="two.sided", var.equal=TRUE) #p=.7172, NS
t.test(but.rc.3,but.rc.28, alternative="two.sided", var.equal=TRUE) #p=0.6331, NS

shapiro.test(but.wd.3) #normal
shapiro.test(but.wd.14) #normal
shapiro.test(but.wd.28) #normal

var.test(x=but.wd.3,y=but.wd.14) #unequal variance
var.test(x=but.wd.3,y=but.wd.28) #equal variance

t.test(but.wd.3,but.wd.14, alternative="greater", var.equal=FALSE) #p=0.06
t.test(but.wd.3,but.wd.28, alternative="greater", var.equal=TRUE) #p=0.181
#but cmon, I think these aren't significant because of a statistical power issue


### Propionate now

prop.rc = subset(prop, Treatment=="RC-ABX")
prop.rc.3 = subset(prop.rc, Rec_day_adj=="-3")
prop.rc.3 = prop.rc.3$value

prop.rc.14 = subset(prop.rc, Rec_day_adj=="14")
prop.rc.14 = prop.rc.14$value

prop.rc.28 = subset(prop.rc, Rec_day_adj=="28")
prop.rc.28 = prop.rc.28$value

prop.wd = subset(prop, Treatment=="WD-ABX")
prop.wd.3 = subset(prop.wd, Rec_day_adj=="-3")
prop.wd.3 = prop.wd.3$value

prop.wd.14 = subset(prop.wd, Rec_day_adj=="14")
prop.wd.14 = prop.wd.14$value

prop.wd.28 = subset(prop.wd, Rec_day_adj=="28")
prop.wd.28 = prop.wd.28$value

shapiro.test(prop.rc.3) #normal
shapiro.test(prop.rc.14) #normal
shapiro.test(prop.rc.28) #normal

var.test(x=prop.rc.3,y=prop.rc.14) #equal variance
var.test(x=prop.rc.3,y=prop.rc.28) #equal variance

t.test(prop.rc.3,prop.rc.14, alternative="two.sided", var.equal=TRUE) #p=.9775, NS
t.test(prop.rc.3,prop.rc.28, alternative="two.sided", var.equal=TRUE) #p=0.895, NS

shapiro.test(prop.wd.3) #not normal
shapiro.test(prop.wd.14) #not normal
shapiro.test(prop.wd.28) #normal

var.test(x=prop.wd.3,y=prop.wd.14) #unequal variance
var.test(x=prop.wd.3,y=prop.wd.28) #equal variance

wilcox.test(prop.wd.3,prop.wd.14, alternative="greater", var.equal=FALSE) #p=0.0218
wilcox.test(prop.wd.3,prop.wd.28, alternative="greater", var.equal=TRUE) #p=0.05
#same as butyrate, i'm guessing it's a statistical power issue

