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
dCt<-cSplit(dCt,"Sample", sep= "_", drop=F)

#rename columns appropriately
dCt<-rename(dCt,replace=c("Sample_1"="Ploidy","Sample_2"="Desiccation","Sample_3"="HeatShock","Sample_4"="Sample"))

#calculate normalized expression of target gene Ct relative to actin Ct using: 2^-(delta Ct)
dCt$HSC70<-2^-(dCt$HSC70-dCt$Actin)
dCt$DNMT1<-2^-(dCt$DNMT1-dCt$Actin)
dCt$MBD2<-2^-(dCt$MBD2-dCt$Actin)
dCt$MeCP2<-2^-(dCt$MeCP2-dCt$Actin)
dCt$HIF1A<-2^-(dCt$HIF1A-dCt$Actin)
dCt$HATHaP2<-2^-(dCt$HATHaP2-dCt$Actin)
dCt$HAT<-2^-(dCt$HAT-dCt$Actin)
dCt$HSP90<-2^-(dCt$HSP90-dCt$Actin)


#log transform the data to develop normality in data
dCt$HSC70log<-log(dCt$HSC70)
dCt$DNMT1log<-log(dCt$DNMT1)
dCt$MeCP2log<-log(dCt$MeCP2)
dCt$HIF1Alog<-log(dCt$HIF1A)
dCt$HATHaP2log<-log(dCt$HATHaP2)
dCt$HATlog<-log(dCt$HAT)
dCt$HSP90log<-log(dCt$HSP90)
dCt$MBD2log<-log(dCt$MBD2)

#Run ANOVA's on all log transformed data as well as Tukey's Honestly Significant Difference post hoc test
HSC70<-aov(HSC70log~Pop+Treat+Pop:Treat, data=dCt)
HSC70
TukeyHSD(HSC70)
summary(HSC70)

DNMT1<-aov(DNMT1log~Pop+Treat+Pop:Treat, data=dCt)
DNMT1
TukeyHSD(DNMT1)
summary(DNMT1)

MeCP2<-aov(MeCP2log~Pop+Treat+Pop:Treat, data=dCt)
MeCP2
TukeyHSD(MeCP2)
summary(MeCP2)

HIF1A<-aov(HIF1Alog~Pop+Treat+Pop:Treat, data=dCt)
HIF1A
TukeyHSD(HIF1A)
summary(HIF1A)

HATHaP2<-aov(HATHaP2log~Pop+Treat+Pop:Treat, data=dCt)
HATHaP2
TukeyHSD(HATHaP2)
summary(HATHaP2)

HAT<-aov(HATlog~Pop+Treat+Pop:Treat, data=dCt)
HAT
TukeyHSD(HAT)
summary(HAT)

HSP90<-aov(HSP90log~Pop+Treat+Pop:Treat, data=dCt)
HSP90
TukeyHSD(HSP90)
summary(HSP90)

PGEEP4<-aov(PGEEP4log~Pop+Treat+Pop:Treat, data=dCt)
PGEEP4
TukeyHSD(PGEEP4)
summary(PGEEP4)

MBD23<-aov(MBD2log~Pop+Treat+Pop:Treat, data=dCt)
MBD23
TukeyHSD(MBD23)
summary(MBD23)

#graph all normalized Ct values to produce boxplots to visualize data

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=HSC70,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.3))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  annotate("text",x=c("C","M","T"), y=0.3, label=c("A", "A", "B"), size=10)+
  labs(x="Treatment", y=expression(paste("HSC70 Expression (",Delta,"Ct)")))


ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=DNMT1, fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.005))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  labs(x="Treatment", y=expression(paste("DNMT1 Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=MeCP2,fill=Pop))+theme_bw()+
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
  labs(x="Treatment", y=expression(paste("MeCP2 Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=HIF1A,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,0.005))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  labs(x="Treatment", y=expression(paste("HIF1A Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=HATHaP2,fill=Pop))+theme_bw()+
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
  labs(x="Treatment", y=expression(paste("HATHaP2 Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=HAT,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,2.5))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  annotate("text",x=c(1.25,2.25), y=2.5, label=c("*","*"), size=12)+
  labs(x="Treatment", y=expression(paste("HAT Expression (",Delta,"Ct)")))

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=HSP90,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,1.5))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  annotate("text",x=c(1.25,2.25), y=1.27, label=c("*","*"), size=12)+
  labs(x="Treatment", y=expression(paste("HSP90 Expression (",Delta,"Ct)")))

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

ggplot(data=dCt)+geom_boxplot(aes(x=Treat, y=MBD2,fill=Pop))+theme_bw()+
  scale_fill_grey(start=0.37, end=.9,
                  labels=c("Dabob Bay","Fidalgo Bay","Oyster Bay"))+
  guides(fill=guide_legend(title="Population"))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25),
        legend.position=c(.09,.87),panel.grid.major=element_blank(),
        legend.key=element_rect(fill=NA))+
  ylim(c(0,.3))+scale_x_discrete(labels=c("Control","Mechanical","Temperature"))+
  annotate("text",x=c("C","M","T"), y=.3, label=c("A", "B", "AB"), size=10)+
  labs(x="Treatment", y=expression(paste("MBD23 Expression (",Delta,"Ct)")))




