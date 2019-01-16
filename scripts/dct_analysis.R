========================
  #
  #UNCOMMENT the lines below if you do have the packages already installed
  #
install.packages("ggplot2")
install.packages("plyr")
install.packages("splitstackshape")


=============================
  
  
  #Necessary Packages to manipulate data and plot values. 
require(plyr)
require(ggplot2)
require(splitstackshape)

#Read in  Ct value table
dCt<-read.csv("data/qpcr_ct_values/qpcr_data_consolidated.csv", header=T)

#Split SAMPLE_ID column to create columns for population, treatment, and sample number
dCt<-cSplit(dCt,"SAMPLE_ID", sep= "_", drop=F)

#rename columns appropriately
dCt<-rename(dCt,replace=c("SAMPLE_ID_1"="Pop","SAMPLE_ID_2"="Treat","SAMPLE_ID_3"="Sample"))

#calculate normalized expression of target gene Ct relative to actin Ct using: 2^-(delta Ct)
dCt$CARM<-2^-(dCt$CarmCt-dCt$Actinct)
dCt$TLR<-2^-(dCt$TLR-dCt$Actinct)
dCt$TRAF<-2^-(dCt$TRAFct-dCt$Actinct)
dCt$H2AV<-2^-(dCt$H2AVct-dCt$Actinct)
dCt$PGRP<-2^-(dCt$PGRP-dCt$Actinct)
dCt$HSP70<-2^-(dCt$HSP70Ct-dCt$Actinct)
dCt$BMP2<-2^-(dCt$BMP2-dCt$Actinct)
dCt$GRB2<-2^-(dCt$GRB2-dCt$Actinct)
dCt$PGEEP4<-2^-(dCt$PGEEP4ct-dCt$Actinct)

#log transform the data to develop normality in data
dCt$CARMlog<-log(dCt$CARM)
dCt$TLRlog<-log(dCt$TLR)
dCt$H2AVlog<-log(dCt$H2AV)
dCt$PGRPlog<-log(dCt$PGRP)
dCt$HSP70log<-log(dCt$HSP70)
dCt$BMP2log<-log(dCt$BMP2)
dCt$GRB2log<-log(dCt$GRB2)
dCt$PGEEP4log<-log(dCt$PGEEP4)
dCt$TRAFlog<-log(dCt$TRAF)

#Run ANOVA's on all log transformed data as well as Tukey's Honestly Significant Difference post hoc test
CARM<-aov(CARMlog~Pop+Treat+Pop:Treat, data=dCt)
CARM
TukeyHSD(CARM)
summary(CARM)

TLR<-aov(TLRlog~Pop+Treat+Pop:Treat, data=dCt)
TLR
TukeyHSD(TLR)
summary(TLR)

H2AV<-aov(H2AVlog~Pop+Treat+Pop:Treat, data=dCt)
H2AV
TukeyHSD(H2AV)
summary(H2AV)

PGRP<-aov(PGRPlog~Pop+Treat+Pop:Treat, data=dCt)
PGRP
TukeyHSD(PGRP)
summary(PGRP)

HSP70<-aov(HSP70log~Pop+Treat+Pop:Treat, data=dCt)
HSP70
TukeyHSD(HSP70)
summary(HSP70)

BMP2<-aov(BMP2log~Pop+Treat+Pop:Treat, data=dCt)
BMP2
TukeyHSD(BMP2)
summary(BMP2)

GRB2<-aov(GRB2log~Pop+Treat+Pop:Treat, data=dCt)
GRB2
TukeyHSD(GRB2)
summary(GRB2)

PGEEP4<-aov(PGEEP4log~Pop+Treat+Pop:Treat, data=dCt)
PGEEP4
TukeyHSD(PGEEP4)
summary(PGEEP4)

TRAF3<-aov(TRAFlog~Pop+Treat+Pop:Treat, data=dCt)
TRAF3
TukeyHSD(TRAF3)
summary(TRAF3)

#graph all normalized Ct values to produce boxplots to visualize data

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=CARM,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.3))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  annotate("text",x=c("C","M","T"), y=0.3, label=c("A", "A", "B"), size=10)+
  labs(x="Treatment", y=expression(paste("CARM Expression (",Delta,"Ct)")))


ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=TLR, fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.005))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  labs(x="Treatment", y=expression(paste("TLR Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=H2AV,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,1.9))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  annotate("text",x=c("C","M","T"), y=1.9, label=c("A", "A", "B"), size=10)+
  annotate("text",x=c(1.25,3.25), y=1.85, label=c("*","*"), size=12)+
  labs(x="Treatment", y=expression(paste("H2AV Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=PGRP,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.005))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  labs(x="Treatment", y=expression(paste("PGRP Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=HSP70,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.80),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,.5))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  annotate("text",x=c("C","M","T"), y=.5, label=c("AB", "A", "B"), size=10)+
  annotate("text",x=c(2.25,3.25), y=.4, label=c("*","*"), size=12)+
  labs(x="Treatment", y=expression(paste("HSP70 Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=BMP2,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,2.5))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  annotate("text",x=c(1.25,2.25), y=2.5, label=c("*","*"), size=12)+
  labs(x="Treatment", y=expression(paste("BMP2 Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=GRB2,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,1.5))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  annotate("text",x=c(1.25,2.25), y=1.27, label=c("*","*"), size=12)+
  labs(x="Treatment", y=expression(paste("GRB2 Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=PGEEP4,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.15))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  labs(x="Treatment", y=expression(paste("PGEEP4 Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=TRAF,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,.3))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  annotate("text",x=c("C","M","T"), y=.3, label=c("A", "B", "AB"), size=10)+
  labs(x="Treatment", y=expression(paste("TRAF3 Expression (",Delta,"Ct)")))




