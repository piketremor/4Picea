#install.packages("devtools")
#install.packages("G:/My Drive/MEForLab",
                 #repos = NULL,
                 #type = "source")
#devtools::install_github("piketremor/MEForLab")

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

setwd("G:/Shared drives/4Picea/4Picea/raw")
setwd("~/Google Drive/Shared drives/4PICEA/4PICEA/raw")
picea <- read.csv("4Picea.csv")
stemform <- read.csv("StemForm.csv")
site <- read.csv("4Picea_30m.csv")

#-------------------------------------------------------------------------------
#join stemform to picea by BLOCK, PLOT, TREE #join site to picea by BLOCK, PLOT
#-------------------------------------------------------------------------------
picea <- picea %>%
  mutate(TREE = as.character(TREE)) %>%
  left_join(stemform %>% mutate(TREE = as.character(TREE)), by = c("BLOCK", "PLOT", "TREE", "SPP")) %>%
  left_join(site, by = c("BLOCK", "PLOT"))

#-------------------------------------------------------------------------------
#clean data
#-------------------------------------------------------------------------------
picea<- picea %>%
  filter(STATUS.23 == 1)
picea$DBH.23 <- as.numeric(picea$DBH.23)
picea$DBH.23[is.na(picea$DBH.23)] <- 0
picea$HT.23 <- as.numeric(picea$HT.23)
picea$uid <- paste0(picea$BLOCK,".",picea$PLOT)

#-------------------------------------------------------------------------------
#calculate the proportion of each species by plot
#-------------------------------------------------------------------------------
species_proportions <- picea %>%
  group_by(BLOCK, PLOT, SPP) %>%
  summarise(Tree_Count = n(), .groups = 'drop') %>%
  group_by(BLOCK, PLOT) %>%
  mutate(Total_Trees = sum(Tree_Count),
         Proportion = Tree_Count / Total_Trees) %>%
  select(BLOCK, PLOT, SPP, Proportion)

print(species_proportions)

picea <- picea %>%
  left_join(species_proportions, by = c("BLOCK", "PLOT", "SPP"))

ggplot(species_proportions, aes(x = factor(PLOT), y = Proportion, fill = SPP)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ BLOCK, scales = "free_x") +
  labs(x = "Plot", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#-------------------------------------------------------------------------------
# calculate bapa and tpa by mixture
#-------------------------------------------------------------------------------
picea$tree.factor <- 10

bapa_tpa_summary <- picea %>%
  filter(STATUS.23 == "1") %>%
  mutate(tree.ba = DBH.23^2 * 0.005454) %>%
  group_by(BLOCK, PLOT, uid, CODE) %>%
  summarise(
    bapa = sum(tree.ba * tree.factor, na.rm = TRUE),  
    tpa = sum(tree.factor, na.rm = TRUE)              
  )

boxplot(bapa~CODE,data=picea)


picea <- picea %>%
  left_join(bapa_tpa_summary, by = c("uid", "CODE", "BLOCK", "PLOT"))

xyplot(bapa~Proportion|CODE,data=picea)

boxplot(bapa~CODE,data=picea)
plot(picea$bapa, picea$tpa)  # most plots are between 20 and 140 bapa

ggplot(picea, aes(x = bapa, y = tpa, color = CODE)) +
  geom_point(alpha = 0.7) +  
  labs(x = "BAPA",
       y = "TPA",
       color = "Species") +
  xlim(20, 140) +  
  theme_minimal()

require(lattice)
xyplot(bapa ~ tpa | CODE, data = picea)
xyplot(bapa ~ tpa | CODE, data = picea, type="l")

#-------------------------------------------------------------------------------
# #DBH distribution
#-------------------------------------------------------------------------------
picea$dbh.class <- 2*as.integer((picea$DBH.23+(2/2))/2)

picea_distribution <- picea %>%
  filter(STATUS.23 == "1") %>%
  mutate(DiameterClass = cut(DBH.23, 
                             breaks = seq(0, 15, by = 1), 
                             right = FALSE)) %>%
  group_by(dbh.class) %>%
  summarise(Count = n(), 
            .groups = 'drop')

ggplot(picea_distribution, aes(x = dbh.class, y = Count, fill = dbh.class)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(title = "Diameter Class Distribution",
       x = "Diameter Classes (DBH in inches)",
       y = "Number of Trees") +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# diameter distribution by CODE
picea$dbh.class <- 2 * as.integer((picea$DBH.23 + (2 / 2)) / 2)

picea_distribution <- picea %>%
  filter(STATUS.23 == "1") %>%
  mutate(DiameterClass = cut(DBH.23, 
                             breaks = seq(0, max(picea$DBH.23, na.rm = TRUE), by = 2), 
                             right = FALSE)) %>%
  group_by(CODE, DiameterClass) %>%
  summarise(Count = n(), 
            .groups = 'drop')

ggplot(picea_distribution, aes(x = DiameterClass, y = Count, fill = DiameterClass)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~ CODE) +  # Create separate plots for each CODE
  labs(x = "Diameter Class (DBH in inches)",
       y = "Number of Trees") +
  scale_x_discrete(labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12"), drop = FALSE) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for better readability

#-------------------------------------------------------------------------------
#generating height model
#-------------------------------------------------------------------------------
require(nlme)

xyplot(HT.23~DBH.23|SPP,data=picea, type="l")
picea$SPP <- as.factor(picea$SPP)
picea$CODE <- as.factor(picea$CODE)
spruce  <- dplyr::filter(picea,SPP=="NS"|SPP=="RS"|SPP=="WS"|SPP=="BS")
xyplot(HT.23~DBH.23|SPP,data=spruce)
str(spruce)

# ht.mod preformed best with AIC of 7029
ht.mod <-  nlme(HT.23 ~ 4.5+exp(a+b/(DBH.23+1)),
                data = spruce,
                fixed = a + b ~ 1,
                #random = a + b ~ 1 | SPP,  # Random intercept and slope for both
                random = a + b ~ 1 | BLOCK/PLOT/SPP,  
                na.action = na.pass,
                start = c(a = 4.5, b = -6),
                control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))
performance(ht.mod)
summary(ht.mod)
ranef(ht.mod)
AIC(ht.mod) #7029

#ht.mod2 <- nlme(HT.23 ~ 4.5+exp(a+b/(DBH.23+1)),
                #data = picea,
                #fixed = a + b ~ 1,
                #random = a + b ~ 1 | BLOCK/PLOT,  # Random intercept and slope for both
                #na.action = na.pass,
                #start = c(a = 4.5, b = -6),
                #control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000)) #AIC 9326

#ht.mod3 <- nlme(HT.23 ~ 4.5+exp(a+b/(DBH.23+1)),
                #data = spruce,
                #fixed = a + b ~ SPP,
                #random = a + b ~ 1 | BLOCK/PLOT,  # Random intercept and slope for both
                #na.action = na.pass,
                #start = c(a = 4.5,0,0,0, b = -6,0,0,0),
                #control = nlmeControl(returnObject = TRUE, msMaxIter = 1000000, maxIter = 500000)) #AIC 7053
#control=lmeControl(opt="optim")


#ht.mod4 <- nlme(HT.23 ~ 4.5+exp((a+b/DBH.23+1)),
                #data = spruce,
                #fixed = a + b ~ SPP,
                #random = a + b ~ 1 | PLOT,  # Random intercept and slope for both
                #na.action = na.pass,
                #start = c(a = 4.5,0,0,0, b = -6,0,0,0),
                #control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000)) #AIC 7350

#ht.mod5 <- nlme(HT.23 ~ 4.5+exp((a+b/DBH.23+1)),
                #data = spruce,
                #fixed = a + b ~ CODE,
                #random = a + b ~ 1 | BLOCK/PLOT,  # Random intercept and slope for both
                #na.action = na.pass,
                #start = c(a = 4.5,0,0,0,0,0,0,0,0,0,0, b = -6,0,0,0,0,0,0,0,0,0,0),
                #control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000)) #AIC 7310

#ht.mod6 <- nlme(HT.23 ~ a + b * DBH.23 * Proportion,
                #data = picea,
                #fixed = a + b ~ 1,
                #random = a + b ~ 1 | PLOT/SPP,
                #na.action = na.pass,
                #start = c(a = 4.5, b = -6),
                #control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000)) #AIC 9697
#AIC(ht.mod, ht.mod2,ht.mod3,ht.mod4,ht.mod5, ht.mod6)
#spruce$ht.fit <- predict(ht.mod2,spruce)
#xyplot(ht.fit~DBH.23|CODE,data=spruce)

# preditc ht.mod on picea to fit all heights, then fit final heights
picea$ht.fit <- predict(ht.mod, picea)
xyplot(ht.fit~DBH.23|SPP,data=picea)
picea$hd <- picea$ht.fit/picea$DBH.23 #height diameter ratio
xyplot(hd~DBH.23|SPP,data=picea)
picea$final.ht <- ifelse(!is.na(picea$HT.23), picea$HT.23, picea$ht.fit)
xyplot(final.ht~DBH.23|SPP,data=picea)
xyplot(final.ht ~ DBH.23 | SPP,
       data = picea,
       subset = SPP %in% c("WS", "BS", "RS", "NS"))

# dealing with outliers in the ht.mod predictions by SPP
# ended up just manually removing the ht from 2 outlier RS stems
#picea$hd2 <- ifelse(picea$hd<3,0,picea$hd)
#picea$hd3 <- ifelse(picea$hd>8,0,picea$hd)
#picea$diff <- picea$hd2-picea$hd3
#picea$hd.new <- 

#-------------------------------------------------------------------------------
#imputing tree heights for all SPP (Wykoff for SPRUCE & FVS for all other SPP)
#-------------------------------------------------------------------------------
#require(MEForLab)

#d <- (ht.mod$coefficients$random)
#ds <- d$PLOT
#BLOCK <- c("B2","B2","B2","B2","B2",
           #"B2","B2","B2","B2","B2","B2",
           #"B3","B3","B3","B3","B3","B3",
           #"B3","B3","B3","B3",
           #"B4","B4","B4","B4","B4",
           #"B4","B4","B4","B4","B4")
#PLOT <- c("1","2","3","4","5","6","7","8","9","10","11",
          #"1","2","3","4","5","6","7","8","9","10",
          #"1","2","3","4","5","6","7","8","9","10")
#re.df <- as.data.frame(ds[,1:2])
#red <- cbind(re.df,BLOCK,PLOT)
#red$uid <- paste0(red$BLOCK,".",red$PLOT)
#red <- red[c(1,2,5)]
#spruces <- left_join(spruce,red,by=c("uid"))
#pairs(spruces[c(28,29,12,23:24,26)])
#plot(spruces$Proportion,spruces$b)
#prop.mod <- lm(a~Proportion,data=spruces)
#summary(prop.mod)

#picea$sp.dum <- ifelse(picea$SPP=="WS"|picea$SPP=="BS"|picea$SPP=="NS"|picea$SPP=="RS",1,0)

# if spruce, use the local model, otherwise use FVS base equations
#picea$local.ht <- ifelse(picea$sp.dum=="1",predict(ht.mod2,picea),0)
#picea$hd <- picea$local.ht/picea$DBH.23
#xyplot(hd~DBH.23|SPP,data=picea,ylim=c(0,20))
#picea$local.ht <- ifelse(picea$hd>12,0,picea$local.ht)
#xyplot(local.ht~DBH.23|SPP,data=picea)
#picea$fvs.ht <- ifelse(picea$local.ht<1,mapply(wykoff.ht,
                                               #SPP=picea$SPP,
                                               #DBH=picea$DBH.23),
                       #0)
#xyplot(fvs.ht~DBH.23|SPP,data=picea)
#picea$HT.23[is.na(picea$HT.23)] <- 0
#picea$model.ht <- ifelse(picea$local.ht>0,picea$local.ht,picea$fvs.ht)
#picea$final.ht <- ifelse(picea$HT.23>0,picea$HT.23,picea$model.ht)

#xyplot(final.ht~DBH.23|SPP,data=picea)
#xyplot(final.ht~DBH.23|CODE,data=picea)


#-------------------------------------------------------------------------------
#basal area larger (bal) - (tree level)
#-------------------------------------------------------------------------------
spruce.bal <- picea%>%
  mutate(basal.area = picea$DBH.23^2*0.005454)%>%
  mutate(bapa = basal.area*10)%>%
  group_by(uid)%>%
  arrange(desc(DBH.23),.by_group = TRUE)%>%
  mutate(bal = lag(cumsum(bapa)))
spruce.bal$bal[is.na(spruce.bal$bal)] <- 0
xyplot(bal~DBH.23|SPP,data=spruce.bal)
xyplot(bal~DBH.23|CODE,data=spruce.bal)

xyplot(bal ~ DBH.23 | CODE,
       data = spruce.bal,
       subset = SPP %in% c("R", "RW", "W", "N","NR","NW","B","BN","BR","BW"))
xyplot(bal ~ DBH.23 | CODE, 
       data = spruce.bal, 
       subset = CODE %in% c("R", "RW", "W", "N", "NR", "NW", "B", "BN", "BR", "BW"),
       layout = c(2, ceiling(length(unique(spruce.bal$CODE))/2))) 

picea <- picea %>%
  left_join(spruce.bal %>% select(uid, BLOCK, PLOT, TREE, bal), 
            by = c("uid", "BLOCK", "PLOT", "TREE"))

#-------------------------------------------------------------------------------
#height larger (htl) - (tree level)
#-------------------------------------------------------------------------------
spruce.htl <- picea%>%
  group_by(uid)%>%
  arrange(desc(final.ht),.by_group = TRUE)%>%
  mutate(htl = lag(cumsum(final.ht)))
spruce.htl$htl[is.na(spruce.htl$htl)] <- 0
xyplot(htl~final.ht|SPP,data=spruce.htl)
xyplot(htl~final.ht|CODE,data=spruce.htl)

xyplot(htl ~ final.ht | SPP,
       data = spruce.htl,
       subset = SPP %in% c("BS", "RS", "NS", "WS"))


picea <- picea %>%
  left_join(spruce.htl %>% select(uid, BLOCK, PLOT, TREE, htl), 
            by = c("uid", "BLOCK", "PLOT", "TREE"))
#-------------------------------------------------------------------------------
#Site Index - (calculated at tree level, summarize by Block and Plot)
#-------------------------------------------------------------------------------
sum(is.na(picea$final.ht))  # Number of missing values in final.ht
picea <- picea %>%
  filter(!is.na(final.ht))

#vicary.site
picea <- picea%>%
  mutate(vicary.si=mapply(vicary.site,SPP="RS",ht=final.ht, age=28))

picea %>%
  group_by(BLOCK, PLOT) %>%
  summarise(mean.vicary.si = mean(vicary.si, na.rm = TRUE))

#steinman.site
picea <- picea%>%
  mutate(steinman.si=mapply(steinman.site, SPP="RS", HT=final.ht, AGE=28))

picea %>%
  group_by(BLOCK, PLOT) %>%
  summarise(mean.steinman.si = mean(steinman.si, na.rm = TRUE))
#-------------------------------------------------------------------------------
#Top Height (vicary height) - (plot level)
#-------------------------------------------------------------------------------
picea <- picea%>%
  mutate(topht=mapply(vicary.height,SPP="RS", age=28, si=vicary.si))

picea %>%
  group_by(BLOCK, PLOT) %>%
  summarise(spruce.topht = mean(topht, na.rm = TRUE))

#-------------------------------------------------------------------------------
#qmd - (plot level)
#-------------------------------------------------------------------------------
spruce.qmd <- picea %>%
  group_by(BLOCK, PLOT) %>%
  summarise(qmd = sqrt(sum(DBH.23^2, na.rm = TRUE) / n()), .groups = "drop")

mean(spruce.qmd$qmd)

picea <- picea %>%
  left_join(spruce.qmd %>% select(BLOCK, PLOT, qmd), 
            by = c("BLOCK", "PLOT"))

#-------------------------------------------------------------------------------
#relative density index - (plot level)
#-------------------------------------------------------------------------------
spruce.rd <- picea %>%
  mutate(rdi = mapply(relative.density.index, bapa = bapa, qmd = qmd))

spruce.rd_summary <- spruce.rd %>%
  group_by(BLOCK, PLOT) %>%
  summarise(rdi = mean(rdi, na.rm = TRUE))  # Ensure only one value per BLOCK-PLOT

picea <- picea %>%
  left_join(spruce.rd_summary, by = c("BLOCK", "PLOT"))

#-------------------------------------------------------------------------------
#relative spacing - (plot level)
#-------------------------------------------------------------------------------
str(relative.spacing)

picea <- picea %>%
  mutate(rs=mapply(relative.spacing,tpa=tpa, domht=topht))

picea %>%
  group_by(BLOCK, PLOT) %>%
  summarise_at(vars(rs), list(name = mean))
#-------------------------------------------------------------------------------
#stand density index - (plot level)
#-------------------------------------------------------------------------------
picea <- picea %>%
  mutate(sdi=mapply(stand.density.index,tpa=tpa, qmd=qmd))

picea %>%
  group_by(BLOCK, PLOT) %>%
  summarise_at(vars(sdi), list(name = mean))

xyplot(qmd~tpa|CODE,data=picea)

#-------------------------------------------------------------------------------
#volume calculation at the individual tree level (ft3)
#-------------------------------------------------------------------------------
picea <- picea%>%
  mutate(vol=mapply(vol.calc,SPP=SPP,DBH=DBH.23,HT=final.ht))
xyplot(vol~DBH.23|CODE,data=picea)

p.vol.code <- ggplot(picea, aes(factor(CODE), vol)) +
  geom_boxplot() +
  ylab("Volume (cubic feet)") +
  xlab("Species Mixture") +
  ggtitle("Volume by Species Mixture")
print(p.vol.code)

# filter for the species of interest
picea.ns <- picea %>%
  filter(SPP == "NS")

# calculate the volume for just NS across plots
picea.ns <- picea.ns %>%
  mutate(vol = mapply(vol.calc, SPP = SPP, DBH = DBH.23, HT = final.ht))

picea.ns<-filter(picea.ns,CODE != "BW"&CODE !="RW"&CODE !="W")
p.vol <- ggplot(picea.ns, aes(factor(CODE), vol)) +
  geom_boxplot()+
  ylab("Volume (cubic feet)") +
  xlab("Species")

p.vol

#-------------------------------------------------------------------------------
# Variable selection procedures using VSURF package for integer response variable
#-------------------------------------------------------------------------------
library(randomForest)
library(VSURF)
library(caTools)
library(pscl)

# premer start 

picea$deductclass <- 5*as.integer((picea$T_DEDUCT+(5/2))/5)
# are we sure we are only fitting trees that have a measurement

hist(picea$deductclass)
picea2 <- dplyr::filter(picea,SPP=="WS"|SPP=="NS"|SPP=="RS"|SPP=="BS")
d.set <- picea2
d.set$T_DEDUCT[is.na(d.set$T_DEDUCT)] <- 999
d.set <- dplyr::filter(d.set,T_DEDUCT<998)
#d.set$HT.23[is.na(d.set$HT.23)] <- 0
#d.set <- dplyr::filter(d.set,HT.23>0)
plot(d.set$DBH.23,d.set$HT.23)

names(d.set)
obs <- d.set$deductclass
#obs[is.na(obs)] <- 0
preds <- d.set[c(4,11,12,21:46,51:54,57:61,65,66,71,76:79)]
#preds[is.na(preds)] <- 0
#vs <- VSURF(preds,obs,ncores = 4)
#vs$varselect.pred
xyplot(deductclass~DBH.23|SPP,data=d.set)
preds2 <- d.set[c(11,12,21:46,51:55,57:61,65,66,71,76:79)]
#preds2[is.na(preds2)] <- 0
#obs <- d.set$deductclass
#obs[is.na(obs)] <- 0
nvs <- VSURF(preds2,obs)
nvs$varselect.pred
names(preds2)

#lai, code, rs, bal, sdi

mod1 <- zeroinfl(deductclass~rs+sdi+bal+roughness+CODE|1,
                   data=d.set,dist="negbin")
summary(mod1)

mod2 <- zeroinfl(deductclass~rs+sdi+bal+roughness+CODE|rs+sdi+bal+roughness+CODE,
                 data=d.set,dist="negbin")
summary(mod2)
AIC(mod1,mod2)

mod3 <- zeroinfl(deductclass~rs+sdi+bal+roughness+CODE|sdi+bal+CODE,
                 data=d.set,dist="negbin")
summary(mod3)
AIC(mod3,mod2)

d.set$wsi <- (d.set$MeanWD-d.set$SWC2)*-1

mod3.1 <- zeroinfl(deductclass~rs+sdi+bal+roughness+CODE|sdi+bal+CODE+SPP,
                 data=d.set,dist="negbin")
summary(mod3.1)
AIC(mod3.1,mod3)

performance(mod3.1)
# dang, that is good. 

picea$fit.deduction <- ifelse(picea$SPP=="RS"|
                              picea$SPP=="WS"|
                              picea$SPP=="BS"|
                              picea$SPP=="NS",predict(mod3.1,type="response"),0)
picea$fit.deduct.class <- (5*as.integer((picea$fit.deduction+(5/2))/5))/100
# will use the modeled reduction
picea$fit.deduct.class[picea$fit.deduct.class>1] <- 1

picea$adj.vol <- picea$vol*picea$fit.deduct.class
plot(picea$vol,picea$adj.vol,ylim=c(0,15),xlim=c(0,15))
picea$final.vol <- picea$vol-picea$adj.vol

#calclate plot level vol 
picea <- picea %>%
  group_by(BLOCK, PLOT) %>%
  mutate(plot.vol = sum(final.vol, na.rm = TRUE))

picea.2 <- picea %>% filter(CODE != "C")
p.vol.code <- ggplot(picea.2, aes(x = factor(CODE), y = plot.vol)) +
  geom_boxplot() +
  ylab("Volume (ft³/ac)") +
  xlab("Species Mixture") +
  theme_minimal(base_size = 14) 
print(p.vol.code)

# for fun, looking at LAI
plot(d.set$wsi,d.set$LAI)
obs.d <- d.set$LAI
preds.d <- d.set[c(12,21:46,51:55,57:61,65,66,71,76:79)]
fun <- VSURF(preds.d,obs.d)

fun$varselect.pred
names(preds.d)

# covariates with LAI - CODE, tri, roughness, qmd, rdi, and nit
d.set$wsi <- (d.set$MeanWD-d.set$SWC2)*-1
plot(d.set$SWC2,d.set$LAI)

ch <- lm(LAI~roughness+SWC2,data=d.set)
summary(ch)
plot(ch)
d.set$fit <- predict(ch,d.set)
d.set$resid <- d.set$LAI-d.set$fit
plot(d.set$fit,d.set$resid)
abline(h=0)

plot(d.set$LAI~d.set$SWC2)
plot(d.set$MeanWD,d.set$LAI)
boxplot(LAI~CODE,data=d.set)

# premer, end. 


str(d.set)
require(pscl)
iz <- zeroinfl(deductclass~SPP+bal+steinman.si+qmd,
               dist="negbin",data=d.set)
summary(iz)

performance(iz)
plot(iz)
izz <- zeroinfl(deductclass~SPP+bal+steinman.si+qmd,
               dist="poisson",data=d.set)
AIC(iz,izz)
summary()

plot(deductclass ~ SPP, data = d.set, subset = deductclass > 0,
     log = "y", main = "Count (positive)")

# SPP, DBH, rs, Densic# SPP, DBH, rs, htl

# let's try something.... 
d.set$deductclass[is.na(d.set$ded)] <- 0
d.set$d.dummy <- as.factor(ifelse(d.set$deductclass>0,1,0))
vs2 <- VSURF(preds,d.set$d.dummy)

vs2$varselect.pred
names(preds)
#sp, qmd , rs

require(pscl)
mod8 <- zeroinfl(deductclass ~dew + qmd + rs,
                 dist = "negbin", data = d.set)
mod8
plot(mod8)
require(DHARMa)

mod9 <- zeroinfl(deductclass ~dew + qmd + rs,
                 , data = d.set)
AIC(mod8,mod9)
summary(mod8)
# the VSURF looks great for the zero inflateed, now for the continuous

much <- dplyr::filter(d.set,deductclass>0)

require(leaps)
much$deductclass <- as.integer(much$deductclass)
plot(much$deductclass)
dredge <- regsubsets(deductclass~SPP+DBH.23+HT.23+LAI+CODE+elevation+
                                   tri+ tpi+roughness+slope+aspect+flowdir+tmin+tmean+ 
                                   tmax+ dew+vpdmax+vpdmin+McNab+Bolstad+Profile+Planform+
                                   Winds10+Winds50+SWI+RAD+ppt+Parent+
                                   dep+ex.k+nit+SWC2+MeanWD+Proportion+Min_depth+WD2000+
                                   WD2020+WHC+ex.mg+ex.ca+ph+bapa+tpa+hd+bal+htl+ 
                                   vicary.si+steinman.si+topht+qmd+rdi+rs+sdi,
                                   data=much,method="exhaustive",really.big=TRUE)


require(performance)
performance(mod8)

logLik(mod8)*-2
BIC(mod8)
mod8null <- update(mod8,.~1)
pchisq(2*logLik(mod8)-logLik(mod8null),df=16,lower.tail = FALSE)/16

expected.counts<-m2$fitted.values
chisq.term<-(obs.counts-expected.counts)^2/expected.counts
df0<-data.frame(k=c(1:length(obs.counts)),Ok=obs.counts,
                Ek=expected.counts,ChiSqk=chisq.term)


#spp, dbh, htl


#-------------------------------------------------------------------------------
# Overyielding & Transgressive Overyielding, looking at the plot level
#-------------------------------------------------------------------------------
# Generate plot estimates
plot_estimates <- picea %>%
  mutate(qmd.2 = qmd(bapa,tpa),
         rd.2 = bapa/sqrt(qmd.2)) %>%
  group_by(BLOCK, PLOT, CODE) %>%  # Include CODE here for species code
  summarise(total_vol = sum(final.vol, na.rm = TRUE),  
            BAPA = mean(bapa),
            TPA = mean(tpa),
            QMD = mean(qmd.2),
            RD = mean(rd.2),
            .groups = 'drop')


# calculate overyielding at plot level
oy <- function(data) {
  results <- data.frame(Mixture = character(), Overyielding = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:nrow(data)) {
    mixture <- as.character(data$CODE[i])  # make CODE a character
    
    if (nchar(mixture) == 2) { # only calculate overyielding for mixtures
      species1 <- substr(mixture, 1, 1)  
      species2 <- substr(mixture, 2, 2)  
      
      volume1 <- data$total_vol[data$CODE == species1]
      volume2 <- data$total_vol[data$CODE == species2]
      
      mixture_volume <- data$total_vol[i]
      
      if (length(volume1) > 0 && length(volume2) > 0) {
        average_volume <- mean(c(volume1, volume2), na.rm = TRUE)
        overyielding <- mixture_volume / average_volume
        
        results <- rbind(results, data.frame(Mixture = mixture, Overyielding = overyielding))
      }
    }
  }
  
  return(results)
}

oy_results <- oy(plot_estimates)
avg_oy <- oy_results %>%
  group_by(Mixture) %>%
  summarise(avg_oy = mean(Overyielding, na.rm = TRUE))
print(avg_oy)

ggplot(avg_oy, aes(x = Mixture, y = avg_oy)) +
  geom_bar(stat = "identity", fill = "grey") +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  
  labs(x = "Species Mixture", 
       y = "Overyielding Threshold") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = c(1))  

        
# calculate transgressive overyielding
calculate_transgressive_overyielding <- function(data) {
  results <- data.frame(Mixture = character(), Transgressive_Overyielding = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:nrow(data)) {
    mixture <- as.character(data$CODE[i]) 
    
    if (nchar(mixture) == 2) {
      species1 <- substr(mixture, 1, 1)
      species2 <- substr(mixture, 2, 2)
      
      volume1 <- data$total_vol[data$CODE == species1]  
      volume2 <- data$total_vol[data$CODE == species2]  
      
      mixture_volume <- data$total_vol[i]  
      
      if (length(volume1) > 0 && length(volume2) > 0) {
        max_volume <- max(volume1, volume2, na.rm = TRUE)  # maximum volume of the two monocultures
        transgressive_overyielding <- mixture_volume / max_volume  # ratio to the larger monoculture volume
        
        results <- rbind(results, data.frame(Mixture = mixture, Transgressive_Overyielding = transgressive_overyielding))
      }
    }
  }
  
  return(results)
}

toy_results <- calculate_transgressive_overyielding(plot_estimates)
avg_toy <- toy_results %>%
  group_by(Mixture) %>%
  summarise(avg_transgressive_oy = mean(Transgressive_Overyielding, na.rm = TRUE))
print(avg_toy)

ggplot(avg_toy, aes(x = Mixture, y = avg_transgressive_oy)) +
  geom_bar(stat = "identity", fill = "grey") +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  
  labs(x = "Species Mixture", 
       y = "Threshold") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = c(1))  

#-------------------------------------------------------------------------------
# bar chart of plot volume by mixture separated by each SPP in the mixture
#-------------------------------------------------------------------------------
picea.vol <- picea %>%
  filter(CODE != "C") %>% 
  group_by(CODE) %>%
  summarise(
    mean.vol = mean(plot.vol, na.rm = TRUE),
    sd.vol = sd(plot.vol, na.rm = TRUE)
  )

picea.vol <- picea.vol %>% #adjust the values for CODE "BR" since something weird was happening before, manually calculated
  mutate(
    mean.vol = ifelse(CODE == "BR", 164.9, mean.vol),  # set BR mean to 164.9
    sd.vol = ifelse(CODE == "BR", 137.3, sd.vol)       # set BR sd to 137.3
  )

ggplot(picea.vol, aes(x = CODE, y = mean.vol, fill = CODE)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean.vol - sd.vol, ymax = mean.vol + sd.vol), width = 0.2) +  # Add error bars
  labs(x = "Species Mixture", y = "Average Plot Volume") +
  theme_minimal() +
  theme(legend.position = "none")
#-------------------------------------------------------------------------------
# ht dbh graphics by each SPP in each mixture it occurs in
#-------------------------------------------------------------------------------
library(patchwork)

# BS plot
plot_bs <- ggplot(picea %>% filter(SPP == "BS", CODE %in% c("B", "BR", "BW", "BN")) %>%
                    group_by(CODE, DBH.23) %>%
                    summarise(avg_final.ht = mean(final.ht, na.rm = TRUE)), 
                  aes(x = DBH.23, y = avg_final.ht, color = CODE)) +
  geom_smooth(method = "loess", aes(group = CODE), se = FALSE, linetype = "solid") +
  labs(x = "DBH (in)", y = "Height (ft)", color = "Species Mixture") +  
  theme_minimal() +
  theme(legend.position = "bottom")

# RS plot
plot_rs <- ggplot(picea %>% filter(SPP == "RS", CODE %in% c("R", "BR", "RW", "NR")) %>%
                    group_by(CODE, DBH.23) %>%
                    summarise(avg_final.ht = mean(final.ht, na.rm = TRUE)), 
                  aes(x = DBH.23, y = avg_final.ht, color = CODE)) +
  geom_smooth(method = "loess", aes(group = CODE), se = FALSE, linetype = "solid") +
  labs(x = "DBH (in)", y = "Height (ft)", color = "Species Mixture") +  
  theme_minimal() +
  theme(legend.position = "bottom")

# NS plot
plot_ns <- ggplot(picea %>% filter(SPP == "NS", CODE %in% c("N", "BN", "NW", "NR")) %>%
                    group_by(CODE, DBH.23) %>%
                    summarise(avg_final.ht = mean(final.ht, na.rm = TRUE)), 
                  aes(x = DBH.23, y = avg_final.ht, color = CODE)) +
  geom_smooth(method = "loess", aes(group = CODE), se = FALSE, linetype = "solid") +
  labs(x = "DBH (in)", y = "Height (ft)", color = "Species Mixture") +  
  theme_minimal() +
  theme(legend.position = "bottom")

# WS plot
plot_ws <- ggplot(picea %>% filter(SPP == "WS", CODE %in% c("W", "BW", "NW", "RW")) %>%
                    group_by(CODE, DBH.23) %>%
                    summarise(avg_final.ht = mean(final.ht, na.rm = TRUE)), 
                  aes(x = DBH.23, y = avg_final.ht, color = CODE)) +
  geom_smooth(method = "loess", aes(group = CODE), se = FALSE, linetype = "solid") +
  labs(x = "DBH (in)", y = "Height (ft)", color = "Species Mixture") +  
  theme_minimal() +
  theme(legend.position = "bottom")

plot_bs + plot_rs + plot_ns + plot_ws + plot_layout(ncol = 2)

#-------------------------------------------------------------------------------
# ht dbh graphics by SPP in each mixture 
#-------------------------------------------------------------------------------
picea.3 <- picea %>%
  filter(SPP %in% c("NS", "RS", "WS", "BS"), 
         CODE != "C",
         !(CODE == "NR" & SPP == "BS")) %>%  
  group_by(SPP, CODE, DBH.23) %>%
  summarise(avg.final.ht = mean(final.ht, na.rm = TRUE), .groups = "drop")

plot.1 <- ggplot(picea.3, aes(x = DBH.23, y = avg.final.ht, color = SPP)) +
  geom_smooth(method = "loess", aes(group = SPP), se = FALSE, linetype = "solid") +
  labs(x = "DBH (in)", y = "Height (ft)", color = "Species") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~CODE, ncol = 2)  

plot.1

#-------------------------------------------------------------------------------
# proportion bapa by the proportion of species in each mixture 
#-------------------------------------------------------------------------------
picea.4 <- picea %>%
  filter(SPP %in% c("NS", "RS", "WS", "BS"), 
         !(CODE %in% c("C", "B", "R", "N", "W")),  # Proper exclusion
         !(CODE == "NR" & SPP == "BS")) %>%  
  group_by(SPP, CODE, Proportion) %>%
  summarise(avg.bapa = mean(bapa, na.rm = TRUE), .groups = "drop")


plot.2 <- ggplot(picea.4, aes(x = Proportion, y = avg.bapa, color = SPP)) +
  geom_smooth(method = "loess", aes(group = SPP), se = FALSE, linetype = "solid") +
  labs(x = "Species Proportion", y = "BAPA", color = "Species") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~CODE, ncol = 2)  

plot.2


picea.4 <- picea %>%
  filter(SPP %in% c("NS", "RS", "WS", "BS"), 
         !(CODE %in% c("C", "B", "R", "N", "W")),  
         !(CODE == "NR" & SPP == "BS")) %>%  
  group_by(SPP, CODE, Proportion) %>%
  summarise(avg.bapa = mean(bapa, na.rm = TRUE), .groups = "drop")

plot.2 <- ggplot(picea.4, aes(x = Proportion, y = avg.bapa, color = SPP)) +
  geom_point(alpha = 0.7, size = 2) +  # Use points instead of lines
  labs(x = "Species Proportion", y = "BAPA", color = "Species") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~CODE, ncol = 2)  

plot.2

#-------------------------------------------------------------------------------
# DBH distribution based on monculture vs. mixed Pretzsch and Biber 2016
#-------------------------------------------------------------------------------
spruce <- picea %>%
  mutate(stand_type = case_when(
    CODE %in% c("NW", "NR", "NB", "RB", "RW", "BW") ~ "Mixed",
    CODE %in% c("N", "W", "B", "R") ~ "Monoculture",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(DBH.23) & !is.na(stand_type))  # Remove rows with NA in DBH.23 or stand_type

ggplot(spruce, aes(x = DBH.23, fill = stand_type, color = stand_type)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  labs(x = "DBH (inches)",
       y = "Density") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(legend.title = element_blank())

ggplot(spruce, aes(x = DBH.23, fill = stand_type, color = stand_type)) +
  geom_histogram(aes(y = ..count..), bins = 30, alpha = 0.5, position = "identity") +  # Histogram with count on y-axis
  labs(x = "DBH (inches)",
       y = "Number of Trees") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(legend.title = element_blank())

#-------------------------------------------------------------------------------
# tree-level final.vol distribution based on monculture vs. mixed (Pretzsch and Biber 2016)
#-------------------------------------------------------------------------------
spruce <- picea %>%
  mutate(
    stand_type = case_when(
      CODE %in% c("NW", "NR", "NB", "RB", "RW", "BW") ~ "Mixed",
      CODE %in% c("N", "W", "B", "R") ~ "Monoculture",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(final.vol) & !is.na(stand_type))  # Remove rows with NA in volume or stand_type

ggplot(spruce, aes(x = final.vol, fill = stand_type, color = stand_type)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  labs(x = "Volume (cubic feet)",
       y = "Density") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(legend.title = element_blank())

ggplot(spruce, aes(x = final.vol, fill = stand_type, color = stand_type)) +
  geom_histogram(aes(y = ..count..), bins = 30, alpha = 0.5, position = "identity") +  # Histogram with count on y-axis
  labs(x = "Volume (cubic feet)",
       y = "Number of Trees") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(legend.title = element_blank())


#-------------------------------------------------------------------------------
# plot.vol distribution based on monculture vs. mixed (Pretzsch and Biber 2016)
#-------------------------------------------------------------------------------
spruce <- picea %>%
  mutate(
    stand_type = case_when(
      CODE %in% c("NW", "NR", "NB", "RB", "RW", "BW") ~ "Mixed",
      CODE %in% c("N", "W", "B", "R") ~ "Monoculture",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(plot.vol) & !is.na(stand_type))  # Remove rows with NA in plot.vol or stand_type

ggplot(spruce, aes(x = plot.vol, fill = stand_type, color = stand_type)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  labs(x = "Plot Volume (cubic feet)",
       y = "Density") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(legend.title = element_blank())

ggplot(spruce, aes(x = plot.vol, fill = stand_type, color = stand_type)) +
  geom_histogram(aes(y = ..count..), bins = 30, alpha = 0.5, position = "identity") +  # Histogram with count on y-axis
  labs(x = "Plot Volume (cubic feet)",
       y = "Number of Trees") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(legend.title = element_blank())

#-------------------------------------------------------------------------------
# now to determine if microsite variations in soil and topographic features influence tree and stand development and lead to o
# final.ht ~ DBH.23 +.... Wykoff Equation 
#-------------------------------------------------------------------------------
hist(picea$final.ht)
picea5 <- dplyr::filter(picea,SPP=="WS"|SPP=="NS"|SPP=="RS"|SPP=="BS")
obs <- picea5$final.ht
#obs[is.na(obs)] <- 0
preds <- picea5[c(4,11,12,21:46,51:55,57:61,65,66,71,76:79)]
#preds[is.na(preds)] <- 0
#vs <- VSURF(preds,obs,ncores = 4)
#vs$varselect.pred
xyplot(final.ht~DBH.23|SPP,data=picea5)
preds2 <- picea5[c(11,12,21:46,51:55,57:61,65,66,71,76:79)]
#preds2[is.na(preds2)] <- 0
#obs <- picea5$final.ht
#obs[is.na(obs)] <- 0
#vs <- VSURF(preds,obs)
#vs$varselect.pred
#names(preds2)

#sdi, qmd, tpa for preds
#rs, bal, bapa for preds2

rf.frame <- ht1$coefficients$random$SPP
rf.frame
rf.frame <- as.data.frame(rf.frame)

library(stringr)
rf.frame2 <- rf.frame %>%
  rownames_to_column(var = "combined_column") %>%
  mutate(
    split_values = str_split(combined_column, "/"),  
    BLOCK = sapply(split_values, `[`, 1),  
    PLOT = as.integer(sapply(split_values, `[`, 2)),  
    SPP = sapply(split_values, `[`, 3)  
  ) %>%
  select(-combined_column, -split_values)  
print(rf.frame2)

picea5 <- picea5 %>%
  left_join(rf.frame2 %>% select(BLOCK, PLOT, SPP, a, b), by = c("BLOCK", "PLOT", "SPP"))


hist(picea5$a)
hist(picea5$b)

summary(picea5$a)
summary(picea5$b)

pairs(~ a + b + rs + bapa + bal, data = picea5)


a1 <- lm(a ~ rs+ bal + qmd, data = picea5)
summary(a1)
AIC(a1)
par(mfrow = c(2,2))
plot(a1)
par(mfrow = c(1,1))

b1 <- lm(b ~ bal, data = picea5)
summary(b1)
AIC(b1)
par(mfrow = c(2,2))
plot(b1)
par(mfrow = c(1,1))

picea5$a.pred <- predict(a1, newdata = picea5)
picea5$b.pred <- predict(b1, newdata = picea5)

#picea5 <- na.omit(picea5)
ht1 <- nlme(HT.23 ~ 4.5 + exp(a.pred + b.pred / (DBH.23 + 1)),
                     data = picea5,  
                     fixed = a.pred + b.pred ~ 1,
                     random = a.pred + b.pred ~ 1 | BLOCK/PLOT/SPP,
                     na.action = na.pass,
                     start = c(a.pred = 4.5, b.pred = -6),
                     control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))

summary(ht1)
AIC(ht1)

picea5 <- picea5 %>%
  mutate(
    BLOCK = as.factor(BLOCK),
    PLOT = as.factor(PLOT),
    SPP = as.factor(SPP)
  )
# orginal ht model
ht1 <- nlme(HT.23 ~ 4.5 + exp(a + b / (DBH.23 + 1)),
            data = picea,
            fixed = a + b ~ 1,
            random = a + b ~ 1 | BLOCK/PLOT/SPP,  
            na.action = na.pass,
            start = c(a = 4.5, b = -6),
            control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))

summary(ht1)
AIC(ht1)


#-------------------------------------------------------------------------------
#bapa ~  
#-------------------------------------------------------------------------------
hist(picea$bapa)
picea5 <- dplyr::filter(picea,SPP=="WS"|SPP=="NS"|SPP=="RS"|SPP=="BS")
obs <- picea5$bapa
#obs[is.na(obs)] <- 0
names(picea5)
preds <- picea5[c(4,6,11,12,21:46,50:61,66,70,71:78,79,85)]
preds2 <- picea5[c(4,6,11,12,21:46,50:61,66,70:72)] #just keep code, spp, site vars
#preds[is.na(preds)] <- 0
vs <- VSURF(preds,obs)
vs$varselect.pred
names(preds)

#rdi for preds, circular
#CODE, SPP, Winds50, Winds10, ex.ca for preds2

#-------------------------------------------------------------------------------
# look at total.vol by plot or tree level?
#-------------------------------------------------------------------------------
hist(picea$plot.vol)
picea5 <- dplyr::filter(picea,SPP=="WS"|SPP=="NS"|SPP=="RS"|SPP=="BS")
obs <- picea5$plot.vol
#obs[is.na(obs)] <- 0
names(picea5)
preds <- picea5[c(4,6,11,12,21:46,50:61,66,70,71:78,79)]
#preds[is.na(preds)] <- 0
vs <- VSURF(preds,obs)
vs$varselect.pred
names(preds)

#preds ex.k, rdi, CODE, qmd, tpa, roughness, LAI, SPP, Planform





