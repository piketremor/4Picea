dev.off()
rm(list=ls())
library(MEForLab)
library(devtools)
library(dplyr)
library(mosaic)
library(forcats)
library(tidyverse)
library(nlme)
library(lme4)
library(performance)
library(ggplot2)
library(MASS)
require(GreenTimbMerch)
require(reshape2)
require(car)

# yield/overyield with pure and mixed species plots
# what is the driver of volume/yield and is there a difference in the treatment?

setwd("C:/Users/ashley.lynn.carter/Documents/GitHub/4Picea/code/raw")
setwd("G:/Shared drives/4Picea/4Picea/raw")
setwd("~/Documents/GitHub/4Picea/code/raw")
setwd("C:/Users/michael.premer/Documents/GitHub/4Picea/code/raw")
picea <- read.csv("4Picea.csv")
picea <- distinct(picea)
stemform <- read.csv("StemForm.csv")
stemform <- distinct(stemform)
site <- read.csv("4Picea_30m.csv")
site <- distinct(site)
site$uid <- paste0(site$BLOCK,".",site$PLOT)
head(picea)
names(picea)

picea<- picea %>%
  filter(STATUS.23 == 1)
picea$DBH.23 <- as.numeric(picea$DBH.23)
picea$DBH.23[is.na(picea$DBH.23)] <- 0
picea$HT.23 <- as.numeric(picea$HT.23)
picea$uid <- paste0(picea$BLOCK,".",picea$PLOT)
picea <- dplyr::filter(picea,SPP=="NS"|SPP=="RS"|SPP=="WS"|SPP=="BS")

# describe the dataset
# table of range of sizes of trees of species within pure and mixed plots

# table of plot and treatment level averages (WSI, LAI, tri, etc.)


# tree level metrics first
ht.mod <-  nlme(HT.23 ~ 4.5+exp(a+b/(DBH.23+1)),
                data = picea,
                fixed = a + b ~ 1,
                #random = a + b ~ 1 | SPP,  # Random intercept and slope for both
                random = a + b ~ 1 | SPP,  
                na.action = na.pass,
                start = c(a = 4.5, b = -6),
                control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))
picea$HT.23[is.na(picea$HT.23)] <- 0
picea$height <- ifelse(picea$HT.23<1,predict(ht.mod,picea,levels=0),picea$HT.23)


picea  <- picea%>%
  filter(.,DBH.23>0)%>%
  mutate(CW = mapply(MCW,SPP=SPP,DBH=DBH.23),
         CA = (CW/2)^2*pi,
         ht.m=height/3.28,
         dbh.cm=DBH.23*2.54,
         ba = DBH.23^2*0.005454,
         vol.m3 = mapply(KozakTreeVol,
                         Bark="ib",
                         SPP=SPP,
                         DBH=dbh.cm,
                         HT=ht.m,
                         Planted=TRUE),
         ccf=CA/43560*100)

tf <- picea
tf$HCB.23[is.na(tf$HCB.23)] <- 0
tfn <- dplyr::filter(tf,HCB.23>0)
tfn$hdr <- tfn$HT.23/tfn$DBH.23

hcb1 <- lm(HCB.23~DBH.23+HT.23+hdr+log(ccf+.1)+SPP,
           data=tfn,na.action=na.exclude)

summary(hcb1)

tfn$X <- predict(hcb1,tfn)

hcb2 <- nls(HCB.23~HT.23/(1+1*exp(-a*X))^1/6,
            data=tfn,
            start=list(a=1))

summary(hcb2)

tfn$lcr <- (tfn$HT.23-tfn$HCB.23)/tfn$HT.23

hcr.mod <- nls(lcr~10*(a/(1+b*ccf)),
               data=tfn,
               start=list(a=1,b=0.005))
summary(hcr.mod)
plot(hcr.mod,ylim=c(-5,5))
performance(hcr.mod)

picea$fit.lcr <- predict(hcr.mod,picea)
picea$fit.HCB <- picea$height-(picea$fit.lcr*picea$height)
picea$HCB.23[is.na(picea$HCB.23)] <- 0
picea$HCB <- ifelse(picea$HCB.23<1,picea$fit.HCB,picea$HCB.23)

plot.sums <- picea%>%
  mutate(ef.ac = 10)%>%
  group_by(BLOCK,PLOT,uid,CODE)%>%
  summarize(bapa = sum(ba*ef.ac),
            tpa = sum(ef.ac),
            qmd = qmd(bapa,tpa),
            CCF = sum(ccf*ef.ac),
            LAI = mean(LAI),
            rd = relative.density.index(bapa,qmd),
            vol.ha = sum(ef.ac*vol.m3)*2.47)

ht40.set <- picea%>%
  group_by(BLOCK,PLOT)%>%
  slice_max(n=4,order_by=DBH.23)%>%
  summarize(ht.40 = mean(height))


ht40.set$SI50 <- mapply(steinman.site,
                        SPP="RS",HT=ht40.set$ht.40,
                        AGE=27)
plot.sums <- left_join(plot.sums,ht40.set)
plot.sums$rs <- relative.spacing(plot.sums$tpa,plot.sums$ht.40)*100

# now do for species, then you are good to join all back to site
sp.sums <- picea%>%
  mutate(ef.ac = 10)%>%
  group_by(BLOCK,PLOT,SPP,CODE)%>%
  summarize(sp.bapa = sum(ba*ef.ac),
            sp.tpa = sum(ef.ac),
            sp.qmd = qmd(sp.bapa,sp.tpa),
            #sp.CCF = sum(ccf*ef.ac),
            sp.rd = relative.density.index(sp.bapa,sp.qmd),
            sp.vol.ha = sum(ef.ac*vol.m3)*2.47)
allin <- left_join(sp.sums,plot.sums)

gg <- allin%>%
  mutate(uid=paste0(BLOCK,".",PLOT))%>%
  group_by(BLOCK,PLOT,SPP,uid)%>%
  summarize(rel.rd = sp.rd/rd,
            rel.vol = sp.vol.ha/vol.ha)


molten <- melt(as.data.frame(gg),id=c("uid","rel.rd","SPP"))
rd.wide <- dcast(molten,uid~SPP,value.var = "rel.rd",mean)
rd.wide <- rename(rd.wide,BS.rd=BS,NS.rd=NS,RS.rd=RS,WS.rd=WS)

fire <- melt(as.data.frame(gg),id=c("uid","rel.vol","SPP"))
vl.wide <- dcast(fire,uid~SPP,value.var="rel.vol",mean)
vl.wide <- rename(vl.wide,BS.rv=BS,NS.rv=NS,RS.rv=RS,WS.rv=WS)

rel.df <- left_join(rd.wide,vl.wide)
fdf <- left_join(plot.sums,site)
fd <- left_join(fdf,rel.df)
head(fd)
fd$uid
fd[is.na(fd)] <- NA


b <- boxcox(lm(vol.ha ~ 1,data=fd))
# Exact lambda
lambda <- b$x[which.max(b$y)]
print(lambda)

fd <- dplyr::filter(fd,CODE!="C")
fd[is.na(fd)] <- 0

modv <- regsubsets(vol.ha~BS_Suitability+WS_Suitability+RS_Suitability+
                     tri+tpi+roughness+SWI+LAI+CCF+rd+
                     slope+aspect+flowdir+RAD+Winds10+SWI+tmean+ppt+
                     WD2000+SWC2+Winds10+Winds50+MeanWD+nit+ex.k+dep+ph+BS.rd+RS.rd+NS.rd+WS.rd+
                     BS.rv+RS.rv+NS.rv+WS.rv,
                   data=fd)
summary(modv)
mod1 <- lm(vol.ha~rd+LAI+dep,data=fd,
           na.action="na.pass")
summary(mod1)

plot(vol.ha~dep,data=fd)
plot(vol.ha~NS.rd,data=fd)
mod2 <- lme(vol.ha~rd+NS.rd+dep,data=fd,
            na.action="na.pass",
            random=~1|BLOCK)

summary(mod2)
plot(mod2,ylim=c(-10,10))
performance(mod2)
vif(mod2)
plot(fd$vol.ha~fd$dep)
fd$vol.fit <- predict(mod2,fd,levels=1)
boxplot(vol.ha~CODE,data=fd)
boxplot(vol.fit~CODE,data=fd)
summary(mod2)

performance(mod2)
# now for overyielding/transgressive overyielding

# overyielding is if the mixture is greater than the weighted average of the best......

fd
head(fd)
names(fd)


keeps <- fd[c(1:4,10,57:60,65)]

head(keeps)
# okay, last try with a nonlinear model. 

blocker <- keeps%>%
  group_by(BLOCK,CODE)%>%
  #rename(MIX=CODE)%>%
  summarize(mean.block.vol = mean(vol.fit))

mix <- left_join(blocker,keeps)
keeps$CODE <- as.factor(keeps$CODE)
separate(keeps$CODE, c("Spec1", "Spec2"))

head(keeps)
head(blocker)

d <- dcast(blocker,BLOCK~CODE)

ke <- left_join(keeps,d)

key <- ke%>%
  mutate(BS.adj = BS.rd*B,
         NS.adj = NS.rd*N,
         RS.adj = RS.rd*R,
         WS.adj = WS.rd*W,
         oy = (BS.adj+NS.adj+RS.adj+WS.adj)/vol.fit)
require(nlme)

keys <- key[c(1:4,12:14,16,17,19,25)]
names(keys)
keys2 <- dplyr::filter(keys,CODE=="BN"|CODE=="BR"|CODE=="BW"|CODE=="NR"|CODE=="NW"|CODE=="RW")
unique(keys2$CODE)

t0 <- lm(oy~CODE,data=keys2)
summary(t0)

t1 <- lme(oy~CODE,
          random=~1|BLOCK,
          data=keys2)

summary(t1)

plot(t1)
performance(t1)
# BN is the only one that over yields, though BR almost does at p = 0.06



pure <- dplyr::filter(blocker,CODE=="B"|CODE=="R"|CODE=="W"|CODE=="N")%>%
  rename(pure.vol=mean.block.vol)
pure$CODE <- as.factor(pure$CODE)
toy.j <- key[c(1:4,10)]
toy.f <- dplyr::filter(toy.j,CODE=="BN"|CODE=="BR"|CODE=="BW"|CODE=="NR"|CODE=="NW"|CODE=="RW")
d <- dcast(pure,BLOCK~CODE)

toys <- left_join(toy.j,d)
head(toys)

tb <- toys%>%
  mutate(N.toy = ifelse(CODE=="NW"||CODE=="NR"||CODE=="BN",vol.fit/N,0),
         W.toy = ifelse(CODE=="BW"||CODE=="RW",vol.fit/W,0),
         B.toy = ifelse(CODE=="BR",vol.fit/B,0),
         toy = (N.toy+W.toy+B.toy))
toy.df <- dplyr::filter(tb,CODE=="BN"|CODE=="BR"|CODE=="BW"|CODE=="NR"|CODE=="NW"|CODE=="RW")

boxplot(toy~CODE,data=toy.df)

toy.1 <- lme(toy~CODE,
             random=~1|BLOCK,
             data=toy.df)
summary(toy.1)

plot(toy.1)



library(ggplot2)
library(ggeffects)
mydf2 <- ggpredict(mod2,terms=c("rd","dep","NS.rd"),
                   type="fixed",
                   ci_level=0.95,
                   interval="confidence",
                   condition=" ")


#png("~/Desktop/SMC_Sinuosity_Model_Ou;tput.png",units='in',height=5.5,width=14,res=1000)
#theme_set(theme_bw(16))

library(ggplot2)
ggplot(mydf2,aes(x=x,y=predicted,colour=group))+
  geom_line(aes(linetype=group,color=group),size=1)+
  labs(x="Relative density",y="Volume (m3) per hectare",panel="Soil depth (cm)")+
  #ylim(0,55)+
  #xlim(4,6)+
  facet_wrap(~facet)+
  theme_bw(18) 
#theme(legend.position="none")
#scale_y_continuous(trans="sq")
#scale_color_manual(values=c('gray0','gray70','gray40'))+
#scale_fill_manual(values=c('gray0','gray70','gray40'), name="fill")
