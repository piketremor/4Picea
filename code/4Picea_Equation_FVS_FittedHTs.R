picea2 <- read.csv("G:/Shared drives/4Picea/4Picea/raw/4Picea.csv")
source('C:/Users/ashley.lynn.carter/Documents/SFR 575 Advanced Biometrics and Modelling/Week 7/furnival.R')


#CleanData
picea2$DBH.23 <- as.numeric(picea2$DBH.23)
picea2$DBH.23[is.na(picea2$DBH.23)] <- 0
picea2$uid <- paste0(picea2$BLOCK,".",picea$PLOT)

library(nlme)

#Chapman-Richards
library(minpack.lm)
m1=nlsLM(HT.23~4.5+b0*(1-exp(-b1*DBH.23))^b2,data=picea2,
         start=c(b0=10.1849,b1=0.00035,b2=0.1122),na.action=na.omit,
         control=nls.lm.control(maxiter = 500))

#Curtis
m2=nlsLM(HT.23~4.5+b0*exp(-b1*DBH.23^b2),data=picea2,
         start=c(b0=10.1849,b1=0.35,b2=-0.1122),na.action=na.omit,
         control=nls.lm.control(maxiter = 500))

#Larsen and Hann
m3=nlsLM(HT.23~4.5+exp(b0+b1*DBH.23^b2),data=picea2,
         start=c(b0=10.1849,b1=0.35,b2=-0.1122),na.action=na.omit,
         control=nls.lm.control(maxiter = 500))
summary(m3)

#Weibull
m4=nlsLM(HT.23~4.5+b0*1-exp(b0+b1*DBH.23^b2),data=picea2,
         start=c(b0=10.1849,b1=0.35,b2=-0.1122),na.action=na.omit,
         control=nls.lm.control(maxiter = 500))

#Wykoff
m5=nlsLM(HT.23~4.5+exp(b0+(b1/DBH.23)),data=picea2,
         start=c(b0=10.1849,b1=0.35),na.action=na.omit,
         control=nls.lm.control(maxiter = 500))

#Power
m6=nlsLM(HT.23~4.5+b0*DBH.23^b1,data=picea2,
         start=c(b0=10.1849,b1=0.35),na.action=na.omit,
         control=nls.lm.control(maxiter = 500))

AIC(m1,m2,m3,m4,m5,m6)
plot(m5)

#Develop mixed-effects model

m1b <- nlme(HT.23~4.5+exp(a+(b/DBH.23+1)),
                data=picea,
                fixed=a+b~1,
                random=a+b~1|BLOCK/PLOT/SPP,
                na.action=na.pass,verbose=T,
                start=c(a=4.5, b=-6),
                control=nlmeControl(returnObject = TRUE,msMaxIter = 10000,maxIter = 5000))
summary(m1b)
furnival(m1b)


#tree.HT$RM=ifelse(tree.HT$FVSspID=='RM',1,0)

m1.lme = lme(I(log(HT.23))~log(DBH.23+1),
             random = ~1|SPP/PLOT/BLOCK,
             data = picea2, 
             weights=varPower(0.2,form=~DBH.23),              
             na.action = na.exclude)
summary(m1.lme)
furnival(m1.lme)
plot(m1.lme)
ranef(m1.lme)

picea2$pHT=predict(m1.lme,newdata=picea2,na.action=na.pass)

#Extract the random effects
library(stringr)
SPP.ranef=ranef(m1.lme)$SPP
PLOT.ranef=ranef(m1.lme)$PLOT
BLOCK.ranef=ranef(m1.lme)$BLOCK

SPP.ranef=data.frame(as.character(rownames(spp.ranef)),spp.ranef)
colnames(SPP.ranef)[1]='SPP'
colnames(SPP.ranef)[2]='b0.SPP'

PLOT.ranef=data.frame(as.character(rownames(mu.ranef)),mu.ranef)
colnames(PLOT.ranef)[1]="PLOT"
colnames(PLOT.ranef)[2]="b0.PLOT"

BLOCK.ranef=data.frame(as.character(rownames(BLOCK.ranef)),BLOCK.ranef)
colnames(BLOCK.ranef)[1]="BLOCK"
colnames(BLOCK.ranef)[2]="b0.BLOCK"

SPP.PLOT=data.frame(str_split_fixed(PLOT.ranef$PLOT, "/", 2))
colnames(SPP.PLOT) <- c("SPP","PLOT")
PLOT.ranef=cbind(PLOT.ranef,SPP.PLOT)
PLOT.ranef=subset(PLOT.ranef,select=c('SPP','PLOT','b0.PLOT'))

SPP.PLOT.BLOCK=data.frame(str_split_fixed(BLOCK.ranef$BLOCK, "/", 3))
colnames(SPP.PLOT.BLOCK) <- c("SPP","PLOT","BLOCK")
BLOCK.ranef=cbind(BLOCK.ranef,SPP.PLOT.BLOCK)
BLOCK.ranef=subset(BLOCK.ranef,select=c('SPP','PLOT','BLOCK','b0.BLOCK'))

#height to crown base
hcb.m1=nlsLM(LLB.HT.23~HT.23/(1+exp(b0+b1*DBH.23^b2)),data=picea2,start=c(b0=0.82,b1=-0.8,b2=-0.49))
hcb.m1a=nlsLM(LLB.HT.23~HT.23/(1+exp(b1*DBH.23^b2)),data=picea2,start=c(b1=-0.8,b2=-0.49))
summary(hcb.m1a)

hcb.m1=nlme(LLB.HT.23~HT.23/(1+exp(b1*DBH.23^b2)),
            data=picea2,fixed=b1+b2~1,
            random=b1~1|SPP/PLOT/BLOCK,
            start=c(b1=0.983,b2=-0.23),
            weights=varPower(0.2,form=~DBH.23),
            na.action=na.omit,
            control=nlmeControl(minScale=1e-6,returnObject=T))
summary(hcb.m1)
furnival(hcb.m1)

library(stringr)
hcb.SPP.ranef=ranef(hcb.m1)$SPP
hcb.PLOT.ranef=ranef(hcb.m1)$PLOT
hcb.BLOCK.ranef=ranef(hcb.m1)$BLOCK

hcb.SPP.ranef=data.frame(as.character(rownames(hcb.SPP.ranef)),hcb.SPP.ranef)
colnames(hcb.SPP.ranef)[1]='SPP'
colnames(hcb.SPP.ranef)[2]='b1.hcb.SPP'

hcb.PLOT.ranef=data.frame(as.character(rownames(hcb.PLOT.ranef)),hcb.PLOT.ranef)
colnames(hcb.PLOT.ranef)[1]="PLOT"
colnames(hcb.PLOT.ranef)[2]="b1.hcb.PLOT"

hcb.SPP.PLOT=data.frame(str_split_fixed(hcb.PLOT.ranef$PLOT, "/", 2))
colnames(hcb.SPP.PLOT) <- c("SPP","PLOT")
hcb.PLOT.ranef=cbind(hcb.PLOT.ranef,hcb.PLOT.mu)
hcb.PLOT.ranef=subset(hcb.PLOT.ranef,select=c('SPP','PLOT','b1.hcb.PLOT'))

##stop

hcb.plt.ranef=data.frame(as.character(rownames(hcb.plt.ranef)),hcb.plt.ranef)
colnames(hcb.plt.ranef)[1]="XXX"
colnames(hcb.plt.ranef)[2]="b1.hcb.plt"

hcb.spp.mu.plt=data.frame(str_split_fixed(hcb.plt.ranef$XXX, "/", 3))
colnames(hcb.spp.mu.plt) <- c("FVSspID","MU","Plot")
hcb.plt.ranef=cbind(hcb.plt.ranef,hcb.spp.mu.plt)
hcb.plt.ranef=subset(hcb.plt.ranef,select=c('FVSspID','MU','Plot','b1.hcb.plt'))

#predict missing values
tree=read.csv('tblUSFSTrees.csv')
tree=tree[!is.na(tree$dbh),]
tree.HT=read.csv('tblUSFSHtCr.csv')
species=read.csv('tblSpCodes.csv')
species$USFSspID=species$PEFspID
tree=merge(tree,species,by=c('USFSspID'))
tree=subset(tree,select=c('MU','Inv','Plot','TreeID','TreeNum','FVSspID','Year','dbh','Expf','Mort','Harv','Ingrowth'))
tree=merge(tree,tree.HT,by=c('MU','Inv','Plot','TreeNum'),all=T)
tree=subset(tree,select=c('MU','Plot','Inv','TreeID','FVSspID','Year','dbh','Expf','Mort','Harv','Ingrowth','HT','BLC'))
tree=merge(tree,spp.ranef,by='FVSspID',all=T)
tree=merge(tree,mu.ranef,by=c('FVSspID','MU'),all=T)
tree=merge(tree,plt.ranef,by=c('FVSspID','MU','Plot'),all=T)
tree$b0.spp=ifelse(is.na(tree$b0.spp),0,tree$b0.spp)
tree$b0.MU=ifelse(is.na(tree$b0.MU),0,tree$b0.MU)
tree$b0.plt=ifelse(is.na(tree$b0.plt),0,tree$b0.plt)

#tree$pHT=4.5+(as.numeric(fixef(m2)[1])+tree$b0.spp+tree$b0.MU)*(1-exp(-as.numeric(fixef(m2)[2])*tree$dbh))^(as.numeric(fixef(m2)[3])+tree$b2.spp+tree$b2.MU)
tree$pHT1=exp(as.numeric(fixef(m1.lme))[1]+tree$b0.spp+tree$b0.MU+tree$b0.plt+as.numeric(fixef(m1.lme))[2]*log(tree$dbh+1))
tree$pHT2=exp(as.numeric(fixef(m1.lme))[1]+tree$b0.spp+tree$b0.MU+tree$b0.plt+as.numeric(fixef(m1.lme))[2]*log(tree$dbh+1)+ + 0.5*sigma(m1.lme)^2)
tree$pHT=ifelse(is.na(tree$HT),1,0)
tree$HT=ifelse(is.na(tree$HT),tree$pHT2,tree$HT)

tree=merge(tree,hcb.spp.ranef,by='FVSspID',all=T)
tree=merge(tree,hcb.mu.ranef,by=c('FVSspID','MU'),all=T)
tree=merge(tree,hcb.plt.ranef,by=c('FVSspID','MU','Plot'),all=T)
tree$b1.hcb.spp=ifelse(is.na(tree$b1.hcb.spp),0,tree$b1.hcb.spp)
tree$b1.hcb.MU=ifelse(is.na(tree$b1.hcb.MU),0,tree$b1.hcb.MU)
tree$b1.hcb.plt=ifelse(is.na(tree$b1.hcb.plt),0,tree$b1.hcb.plt)
tree$HCB.ft=tree$HT/(1+exp((as.numeric(fixef(hcb.m1)[1])+tree$b1.hcb.spp+tree$b1.hcb.MU+tree$b1.hcb.plt)*tree$dbh^as.numeric(fixef(hcb.m1)[2])))
tree$p.HCB=ifelse(is.na(tree$BLC),1,0)
tree$BLC=ifelse(is.na(tree$BLC),tree$HCB.ft,tree$BLC)

#compare to FVS
source('fvsne_eqs.r')
tree2=read.csv('tblUSFSTrees.csv')
tree.HT=read.csv('tblUSFSHtCr.csv')
species=read.csv('tblSpCodes.csv')
species$USFSspID=species$PEFspID

tree2=merge(tree2,species,by=c('USFSspID'))
tree.HT=merge(tree.HT,tree2,by=c('MU','Plot','Inv','TreeNum'))
tree.HT$CR=((tree.HT$HT-tree.HT$BLC)/tree.HT$HT)*100
tree.HT=tree.HT[tree.HT$HT<200 & !is.na(tree.HT$FVSspID),]

aa=sapply(tree.HT$FVSspID,FVSNE.SPP)

tree.HT$HTDBH_EQ=t(aa)[,1]
tree.HT$Dbw=as.numeric(t(aa)[,2])
tree.HT$SPPGR=as.numeric(t(aa)[,3])

tree.HT$HT.FVS=mapply(FVSNE.HTDBH,SPP=tree.HT$FVSspID,HTDBH_EQ=tree.HT$HTDBH_EQ,DBH=tree.HT$dbh,Dbw=tree.HT$Dbw)
tree.HT$CR.FVS=mapply(FVSNE.CR,SPPGR=tree.HT$SPPGR,BA=60,DBH=tree.HT$dbh)
tree.HT$BLC.FVS=(tree.HT$CR.FVS/100)*tree.HT$HT

plot(tree.HT$HT.FVS,tree.HT$HT)
abline(a=0,b=1)

plot(tree.HT$CR.FVS,tree.HT$CR)
abline(a=0,b=1)

plot(tree.HT$BLC.FVS,tree.HT$BLC)
abline(a=0,b=1)

FVS.HT.m1=lme(HT~HT.FVS,data=tree.HT,random=~1|Type/FVSspID/MU/Plot/Inv,
              na.action=na.omit)
summary(FVS.HT.m1)
furnival(FVS.HT.m1)

FVS.BLC.m1=lme(BLC~BLC.FVS,data=tree.HT,random=~1|Type/FVSspID/MU/Plot/Inv,
               na.action=na.omit)
summary(FVS.BLC.m1)
furnival(FVS.BLC.m1)

tree.HT$HT.FVS.cal=predict(FVS.HT.m1,newdata=tree.HT,na.action=na.pass)
tree.HT$BLC.FVS.cal=predict(FVS.BLC.m1,newdata=tree.HT,na.action=na.pass)

plot(tree.HT$HT.FVS.cal,tree.HT$HT)
abline(a=0,b=1)

plot(tree.HT$BLC.FVS.cal,tree.HT$BLC)
abline(a=0,b=1)
