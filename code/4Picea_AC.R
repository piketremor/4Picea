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

#min(site$elevation, na.rm = TRUE)  # Minimum elevation
#mean(site$elevation, na.rm = TRUE) # Mean elevation
#max(site$elevation, na.rm = TRUE)
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
species.prop <- picea %>%
  group_by(BLOCK, PLOT, SPP) %>%
  summarise(tree.ct = n(), .groups = 'drop') %>%
  group_by(BLOCK, PLOT) %>%
  mutate(t.trees = sum(tree.ct),
         Proportion = tree.ct / t.trees) %>%
  select(BLOCK, PLOT, SPP, Proportion)

print(species.prop)

picea <- picea %>%
  left_join(species.prop, by = c("BLOCK", "PLOT", "SPP"))

ggplot(species.prop, aes(x = factor(PLOT), y = Proportion, fill = SPP)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ BLOCK, scales = "free_x") +
  labs(x = "Plot", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#-------------------------------------------------------------------------------
# calculate bapa and tpa by mixture
#-------------------------------------------------------------------------------
picea$tree.factor <- 10

bapa.tpa.summary <- picea %>%
  filter(STATUS.23 == "1") %>%
  mutate(tree.ba = DBH.23^2 * 0.005454) %>%
  group_by(BLOCK, PLOT, uid, CODE) %>%
  summarise(
    bapa = sum(tree.ba * tree.factor, na.rm = TRUE),  
    tpa = sum(tree.factor, na.rm = TRUE)              
  )

picea <- picea %>%
  left_join(bapa.tpa.summary, by = c("uid", "CODE", "BLOCK", "PLOT"))

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
xyplot(bapa~Proportion|CODE,data=picea)

#sums <- picea %>%
  #group_by(CODE) %>%
  #summarise(
    #min = min(sdi),
    #mean = mean(sdi),
    #max = max(sdi),
    #sd = sd(sdi)
  #)

#print(sums, n = nrow(sums))

#sums2 <- picea %>%
  #group_by(CODE, SPP) %>%
  #summarise(
    #min = min(final.hcb),
    #mean = mean(final.hcb),
    #max = max(final.hcb),
    #sd = sd(final.hcb)
  #

#print(sums2, n = nrow(sums))
#print(sums2, n = Inf)

#convert to metric 
#error because need to calculate qmd, vol first which is further down
#picea.sum <- picea %>%
  #mutate(
    #qmdcm = qmd * 2.54,  
    #baph = bapa * 2.47105,  
    #tph = tpa * 2.47105,
   #plotvolha = plot.vol * 0.0693 
  #)

#uniquedata <- picea.sum %>%
 #distinct(bapa, BLOCK, PLOT, CODE)

#print(uniquedata, n = Inf)

#avg <- picea %>%
  #group_by(CODE) %>%
  #summarise(
    #avg_baph = mean(topht, na.rm = TRUE),
    #.groups = 'drop'  
  #)

#print(avg, n = Inf)  

#-------------------------------------------------------------------------------
# #DBH distribution
#-------------------------------------------------------------------------------
picea$dbh.class <- 2*as.integer((picea$DBH.23+(2/2))/2)

picea.dist <- picea %>%
  filter(STATUS.23 == "1") %>%
  mutate(DiameterClass = cut(DBH.23, 
                             breaks = seq(0, 15, by = 1), 
                             right = FALSE)) %>%
  group_by(dbh.class) %>%
  summarise(Count = n(), 
            .groups = 'drop')

ggplot(picea.dist, aes(x = dbh.class, y = Count, fill = dbh.class)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(title = "Diameter Class Distribution",
       x = "Diameter Classes (DBH in inches)",
       y = "Number of Trees") +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# diameter distribution by CODE
picea$dbh.class <- 2 * as.integer((picea$DBH.23 + (2 / 2)) / 2)

picea.dist <- picea %>%
  filter(STATUS.23 == "1") %>%
  mutate(DiameterClass = cut(DBH.23, 
                             breaks = seq(0, max(picea$DBH.23, na.rm = TRUE), by = 2), 
                             right = FALSE)) %>%
  group_by(CODE, DiameterClass) %>%
  summarise(Count = n(), 
            .groups = 'drop')

ggplot(picea.dist, aes(x = DiameterClass, y = Count, fill = DiameterClass)) +
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
#generating crown ratio model
#-------------------------------------------------------------------------------
#library(dplyr)
#library(nlme)

#crownratio  <- dplyr::filter(picea,SPP=="NS"|SPP=="RS"|SPP=="WS"|SPP=="BS")
#xyplot(HCB.23~DBH.23|SPP,data=crownratio)


# Northeast FVS Variant Guide
# WS SPP Group 3, RS BS NS SPP Group 4
# FVS NE Variant Coef by SPP Group
#coefficients <- data.frame(
#sppgroup = c(3, 4),
#b1 = c(7.840, 5.540),
#b2 = c(0.0057, 0.0072),
#b3 = c(1.272, 4.200),
#b4 = c(-0.1420, -0.0530)
#)

#crownratio <- crownratio %>%
#mutate(sppgroup = case_when(
#SPP == "WS" ~ 3,  
#SPP %in% c("NS", "BS", "RS") ~ 4,  
#TRUE ~ NA_integer_  # NA for other SPP
#))

#crownratio <- crownratio %>%
  #left_join(coefficients, by = "sppgroup")  



# impute missing crown ratios using the basic formula from FVS NE Variant
# guessing this doesn't calibrate to the LCR.23 heights measured in the field
# looks like this assigned the same cr.fit value by SPP don't vary that much .2-.5 CR. For example, I had a measured CR on 0.82 but the predicted CR was .28. 
#crownratio <- crownratio %>%
  #mutate(
  #cr.fit = (10 * b1 / (1 + b2 * bapa) + (b3 * (1 - exp(-b4 * DBH.23)))) / 100,  
  #final.cr = ifelse(!is.na(LCR.23), LCR.23, cr.fit)  
  #)

#library(lattice)
#xyplot(final.cr ~ DBH.23 | SPP, data = crownratio, subset = SPP %in% c("WS", "BS", "RS", "NS"))

# crown ratio model using nlme/following ht.mod layout
# could not get to run, is the issue because of the 2 species groups so the start parameters vary?
#cr.mod <- nlme(
  #LCR.23 ~ 10 * (b1 / (1 + b2 * bapa) + (b3 * (1 - exp(-b4 * DBH.23)))),  
  #data = crownratio,  
  #fixed = b1 + b2 + b3 + b4 ~ 1,
  #random = b1 + b2 + b3 + b4 ~ 1 | BLOCK/SPP,  
  ##start = c(b1 = 7.840, b2 = 0.0057, b3 = 1.272, b4 = -0.1420))#, 
  #control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000, optimMethod = "L-BFGS-B"))

#plot(crownratio$DBH.23,crownratio$HCB.23)
#xyplot(HCB.23~DBH.23|SPP,data=crownratio)
#crownratio$CL <- crownratio$HT.23-crownratio$HCB.23
#xyplot(CL~HT.23|SPP,data=crownratio,type=c("p"))
#boxplot(CL~CODE,data=crownratio)

#m2 <- lm(LCR.23~CODE+SPP,data=crownratio)
#summary(m2)


#m3 <- lm(LCR.23~CODE+SPP+DBH.23+HT.23+log(DBH.23),data=crownratio)
#summary(m3)
#m4 <- lme(LCR.23~CODE+SPP+DBH.23+HT.23+log(DBH.23)+MeanWD,data=crownratio,random=~1|BLOCK,na.action="na.omit")

#crownratio$fit <- predict(m3,crownratio)
#xyplot(crownratio~HT.23|SPP,data=crownratio)
#Error in `[.data.frame`(y, id) : undefined columns selected
#xyplot(fit~HT.23|SPP,data=crownratio)
#xyplot(fit~HT.23|SPP,data=crownratio,type="l")
#xyplot(fit~HT.23|SPP,data=crownratio,type="p")
#xyplot(1-fit~HT.23|SPP,data=crownratio,type="p")
#crownratio$crowndepth <- (crownratio$HT.23*crownratio$fit)
#xyplot(crowndepth~HT.23|SPP,data=crownratio,type="p")
#crownratio$crowndepth <- crownratio$HT.23-(crownratio$HT.23-(crownratio$HT.23*crownratio$fit))
#xyplot(crowndepth~HT.23|SPP,data=crownratio,type="p")


#cr.mod2 <-  nlme(HCB.23 ~ 4.5+exp(a+b/(DBH.23+1)),
                #data = crownratio,
                #fixed = a + b ~ 1,
                #random = a + b ~ 1 | BLOCK/PLOT/SPP,  
                #na.action = na.pass,
                #start = c(a = 4.5, b = -6),
                #control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))

#performance(cr.mod2)
#summary(cr.mod2)
#ranef(cr.mod2)
#AIC(cr.mod2)

#picea$cr.fit2 <- predict(cr.mod2, crownratio)
#xyplot(cr.fit2 ~ DBH.23 | SPP, data = crownratio)

#picea$final.cr2 <- ifelse(!is.na(picea$LCR.23), picea$LCR.23, picea$cr.fit2)
#xyplot(final.cr2 ~ DBH.23 | SPP, data = crownratio)
#xyplot(final.cr2 ~ DBH.23 | SPP, data = crownratio, subset = SPP %in% c("WS", "BS", "RS", "NS"))

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
#maximum crown width (tree level)
#-------------------------------------------------------------------------------
picea <- picea %>%
  mutate(DBH = as.numeric(DBH.23), 
         MCW = mapply(MCW, SPP = SPP, DBH = DBH))


#-------------------------------------------------------------------------------
#CCF 
#-------------------------------------------------------------------------------
#picea <- picea %>%
  #mutate(
    #CL = HT.23 - HCB.23,   #CL = crown length
    #CD = (CL + MCW) / 2     # CD = crown width diameter
  #)

#picea <- picea %>%
  #mutate(
    #CA = pi * (CD / 2)^2  # CA = crown area
  #)

#picea <- picea %>%
  #group_by(BLOCK, PLOT) %>%
  #mutate(CA.plot = sum(CA, na.rm = TRUE)) %>%
  #ungroup()

#picea <- picea %>%
  #mutate(CCF = (CA / CA.plot) * 10) #1/10 acre plots


# mikes formula
spruce.ccf <- picea%>%
  mutate(CW = mapply(MCW,SPP=SPP,DBH=DBH),
         CA = (CW/2)^2*pi,
         CA.exp = CA*10)%>%
  group_by(BLOCK,PLOT)%>%
  summarize(CCF = (sum(CA.exp)/43560)*100) #should I divide by 4356 to represent 1/10 acre plots or no

#picea.m <- left_join(spruce,spruce.ccf)
#plot(picea.m$bapa,picea.m$CCF)
  
picea <- picea %>%
  left_join(spruce.ccf %>% select(BLOCK, PLOT, CCF), 
            by = c("BLOCK", "PLOT"))

--------------------------------------------------------------------
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
#Top Height/HT40 (vicary height) - (tree or plot level?)
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

#calculate volume at plot level
picea <- picea %>%
  group_by(BLOCK, PLOT) %>%
  mutate(stand.vol = sum(vol, na.rm = TRUE))

# predict vol from LAI? roughness?
vol <- lm(stand.vol ~ LAI + roughness + factor(CODE), data = picea)
summary(vol)
AIC(vol)

vol2 <- lm(stand.vol ~ I(log(LAI)) + I(log(roughness)) + factor(CODE), data = picea)
summary(vol2)
AIC(vol2)

vol3 <- lm(stand.vol ~ roughness + I(log(LAI)) + factor(CODE), data = picea)
summary(vol3)
AIC(vol3)

#plot(vol)

plot(vol3$fitted.values, resid(vol), 
     xlab = "Fitted Values", 
     ylab = "Residuals", 
     main = "Residual Plot",
     pch = 20, col = "black")
abline(h = 0, lty = 2, col = "red")

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
#nvs <- VSURF(preds2,obs)
#nvs$varselect.pred
#names(preds2)

#lai, code, rs, bal, sdi

#mod1 <- zeroinfl(deductclass~rs+sdi+bal+roughness+CODE|1,
                   #data=d.set,dist="negbin")
#summary(mod1)

#mod2 <- zeroinfl(deductclass~rs+sdi+bal+roughness+CODE|rs+sdi+bal+roughness+CODE,
                 #data=d.set,dist="negbin")
#summary(mod2)
#AIC(mod1,mod2)

#mod3 <- zeroinfl(deductclass~rs+sdi+bal+roughness+CODE|sdi+bal+CODE,
                 #data=d.set,dist="negbin")
#summary(mod3)
#AIC(mod3,mod2)

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

#convert to metric 
#picea.2 <- picea.2 %>%
  #mutate(plot.vol_m3ha = plot.vol * 0.070)

#p.vol.code <- ggplot(picea.2, aes(x = factor(CODE), y = plot.vol_m3ha)) +
  #geom_boxplot() +
  #ylab("Volume (m³/ha)") +
  #xlab("Species Mixture") +
  #theme_minimal(base_size = 14) 

#print(p.vol.code)

# for fun, looking at LAI
plot(d.set$wsi,d.set$LAI)
obs.d <- d.set$LAI
preds.d <- d.set[c(12,21:46,51:55,57:61,65,66,71,76:79)]
#fun <- VSURF(preds.d,obs.d)

#fun$varselect.pred
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
                 data = d.set)
AIC(mod8,mod9)
summary(mod8)
# the VSURF looks great for the zero inflateed, now for the continuous

much <- dplyr::filter(d.set,deductclass>0)

require(leaps)
much$deductclass <- as.integer(much$deductclass)
plot(much$deductclass)
#dredge <- regsubsets(deductclass~SPP+DBH.23+HT.23+LAI+CODE+elevation+
 #                                  tri+ tpi+roughness+slope+aspect+flowdir+tmin+tmean+ 
  #                                 tmax+ dew+vpdmax+vpdmin+McNab+Bolstad+Profile+Planform+
   #                                Winds10+Winds50+SWI+RAD+ppt+Parent+
    #                               dep+ex.k+nit+SWC2+MeanWD+Proportion+Min_depth+WD2000+
     #                              WD2020+WHC+ex.mg+ex.ca+ph+bapa+tpa+hd+bal+htl+ 
      #                             vicary.si+steinman.si+topht+qmd+rdi+rs+sdi,
       #                            data=much,method="exhaustive",really.big=TRUE)


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
# HCB Model
#-------------------------------------------------------------------------------
#h.set <- picea10
#h.set$HCB.23[is.na(h.set$HCB.23)] <- 999
#h.set <- dplyr::filter(h.set,HCB.23<998)

#names(h.set)
#obs <- h.set$HCB.23
#obs[is.na(obs)] <- 0
#preds <- h.set[c(4,6,7,11,12,21:46,50:61,65,66,71,72,75,76,77,78,79,80,81,82)]
#preds2 <- h.set[c(4,6,7,11,12,21:46,50:54,57:61)] #take out stand structure metrics
#preds[is.na(preds)] <- 0
#vs <- VSURF(preds2,obs,ncores = 4)
#vs$varselect.pred
#names(preds2)

#HT.23, log(bal), log(CCF), hd, DBH.23 for A. Weiskittel paper
#SPP, CODE, rdi, vicary.si (circular), bal for preds
#SPP, HT.23, DBH.23, CODE for preds2
# tried just leaving site vars in as only predictors but only picked up RS_Suitability

# your filter below includes trees with no ht measurent, right? 

#picea.11 <- dplyr::filter(picea,SPP=="WS"|SPP=="NS"|SPP=="RS"|SPP=="BS")
#picea.11 <- dplyr::filter(picea, (SPP %in% c("WS", "NS", "RS", "BS")) & (CODE != "C")) # filter out the control (C) plots too?
#names(spruce.bal)
#spruce.bal2 <- spruce.bal[c(1,2,3,72)]

d.set <- dplyr::left_join(d.set,spruce.bal2)
xyplot(HCB.23~HT.23|SPP,data=d.set)

#hcb.mod <- lm(HCB.23 ~I(log(CCF)) + final.ht + bal + factor(SPP) + factor(CODE), data = d.set)
#summary(hcb.mod)

#hcb.mod2 <- lm(HCB.23 ~final.ht + bal  + factor(SPP) + factor(CODE),data = d.set)
#summary(hcb.mod2)

#hcb.mod3 <- lme(HCB.23~final.ht+bal+factor(SPP)+factor(CODE),data=d.set,

#model2 <- lm(HCB.23 ~I(log(CCF))+HT.23 + bal  + factor(SPP) + factor(CODE),data = d.set)
#summary(model2)

#model3 <- lme(HCB.23~HT.23+I(log(CCF))+bal+factor(SPP)+factor(CODE),data=d.set,
          #random=~1|BLOCK/PLOT,na.action="na.omit",method="REML")


hcb.mod4 <- lm(HCB.23 ~final.ht + DBH.23  + factor(SPP) + factor(CODE),data = d.set) #log final.ht or DBH.23 did not improve AIC
summary(hcb.mod4)
AIC(hcb.mod4)
#AIC(model2, hcb.mod4) #hcb.mod4 is better model 

residuals <- residuals(hcb.mod4)
ggplot(data.frame(fitted = fitted(hcb.mod4), residuals = residuals(model)), aes(x = fitted, y = residuals)) +
  geom_point(color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

rmse <- rmse(hcb.mod4)
print(paste("RMSE:", rmse))

MB <- mean(residuals)
print(paste("Mean Bias (MB):", MB))

MAB <- mean(abs(residuals))
print(paste("Mean Absolute Bias (MAB):", MAB))

# predict hcb.mod4 on all spruce, and join to picea df
picea$hcb.fit <- ifelse(picea$SPP=="RS"|
                          picea$SPP=="WS"|
                          picea$SPP=="BS"|
                          picea$SPP=="NS",predict(hcb.mod4,type="response"),0)
picea$final.hcb <- ifelse(!is.na(picea$HCB.23), picea$HCB.23, picea$hcb.fit)


#now, for the model. 

# test
#d.set$fit.X <- predict(hcb.mod4,d.set)

#d.set$fit.hcb <- d.set$HT.23/
  #(((1+1*exp(-1*d.set$fit.X)))^(1/6))

#plot(d.set$fit.hcb,d.set$HCB.23)
#abline(0,1)

#xyplot(fit.hcb~HCB.23|CODE,data=d.set,type="l")
#xyplot(fit.hcb~HT.23|CODE,data=d.set,type="l")

#d.set$lin.fit <- predict(model2,d.set)
#plot(d.set$HCB.23,d.set$lin.fit)
#abline(0,1)

#xyplot(HCB.23~lin.fit|SPP,data=d.set,xlim=c(0,30),ylim=c(0,40))

#plot(residuals(model2))

#d.set$exp.fit <- d.set$HT.23*(1-1*exp(-1*d.set$lin.fit^10))
#plot(d.set$HCB.23,d.set$exp.fit)
#abline(0,1)


#xyplot(d.set$HCB.23~d.set$DBH.23|d.set$SPP)


# model2 is better
#coefficients <- coef(model2)
#print(coefficients)

#B0 <- -0.27150984       # Intercept
#B1 <- 0.33969644     # HT.23
#B2 <- 0.02672539     #bal  
#B3 <- -4.3512    # factor(SPP)NS
#B4 <- -3.26968020  # factor(SPP)RS
#B5 <- -2.48090535     # factor(CODE)BR
#B6 <- 1.55275136    # factor(CODE)BW
#B7 <- 1.67330013    # factor(CODE)N
#B8 <- 2.02087534     # factor(CODE)NW
#B9 <- -1.57597904    #factor(CODE)R

#d.set$X <- B0 + 
  #B1 +#* d.set$HT.23 + 
  #B2 +#* d.set$bal + 
  #ifelse(d.set$SPP == "NS", B3, 0) + 
  #ifelse(d.set$SPP == "RS", B4, 0) +
  # ifelse(d.set$CODE == "BR", B5, 0) +
  #ifelse(d.set$CODE == "BW", B6, 0) +
  #ifelse(d.set$CODE == "N", B7, 0) +
  #ifelse(d.set$CODE == "NW", B8, 0) +
  #ifelse(d.set$CODE == "R", B9, 0)


#c <- 1
#k <- 1
#m <- 10

#d.set <- d.set %>%
  #mutate(HCB1 = HT.23 / (1 + c * exp(-k * X))^(1/m),
         #HCB2 = HT.23*(1-1*exp(-1*d.set$X^10)))

#plot(d.set$HCB1,d.set$HCB.23)
#abline(0,1)

#plot(d.set$HCB.23,d.set$HCB1) 
#abline(0,1)
#d.set$HT.23*(1-1*exp(-1*d.set$lin.fit^10))

picea$fit.hcb <- predict(hcb.mod4,picea)
spruce.only <- dplyr::filter(picea,SPP=="RS"|SPP=="NS"|SPP=="WS"|SPP=="BS")
spruce.only$HT.23 <- ifelse(is.na(spruce.only$HT.23),spruce.only$final.ht,spruce.only$HT.23)
spruce.only$fit.hcb <- predict(hcb.mod4,spruce.only)
spruce.only$fit.hcb <- ifelse(spruce.only$fit.hcb<0.3,0.3,spruce.only$fit.hcb)

cpi.frame <- spruce.only%>%
  group_by(BLOCK,PLOT)%>%
  summarize(Top = max(HT.23),
            Low = min(fit.hcb))
cpi.frame <- as.data.frame(cpi.frame)

cpi.frame$uid.b <- paste0(cpi.frame$BLOCK,".",cpi.frame$PLOT)
nd <- data.frame(uid.b = rep(c(cpi.frame$uid.b),each=100),
                 prop = seq(from=0,to=1,length=100))
nd2 <- left_join(nd,cpi.frame)
head(nd2)

names(spruce.only)
al <- spruce.only[c(1:8,12,66,70,91)]
head(al)
pa <- left_join(nd2,al)
head(pa)
pa$top.prop <- pa$final.ht/pa$Top
pa$bot.prop <- pa$fit.hcb/pa$Top
pa$crown.point <- ifelse(pa$prop>=pa$bot.prop&pa$prop<=pa$top.prop,1,0)
head(pa)
xyplot(crown.point~DBH.23|SPP,data=pa)
#pa$prop <- as.integer(pa$prop)
cr.demog <- pa%>%
  mutate(ef = 10)%>%
  group_by(BLOCK,PLOT,SPP,CODE,prop)%>%
  summarize(cr.points = sum((crown.point)*10))

COLORS = c("blue", "gold", "red", "green")


xyplot(cr.points~prop|CODE,data=cr.demog,type=c("p","l"),
       group=SPP,auto.key=TRUE)
# obviously the spruce only filter doesn't work.. oh well. 
ca <- dplyr::filter(cr.demog,cr.points>0&CODE=="RW"|CODE=="NR"|CODE=="NW"|CODE=="BR"|CODE=="BW"|CODE=="BN")

ca <- dplyr::filter(cr.demog, cr.points > 0 & 
                      (CODE %in% c("RW", "NR", "NW", "BR", "BW", "BN")) & 
                      (SPP %in% c("BS", "NS", "RS", "WS")) &
                      !((SPP %in% c("BS", "WS") & CODE == "NR") |  
                          (SPP %in% c("RS", "WS") & CODE == "BN") |
                          (SPP %in% c("BS", "NS") & CODE == "RW") |
                          (SPP %in% c("RS", "NS") & CODE == "BW")))  

ca$SPP <- factor(ca$SPP, levels = c("BS", "NS", "RS", "WS"))

xyplot(prop~cr.points|CODE,data=ca,type=c("l"))


ca.avg <- ca %>%
  group_by(CODE, SPP, prop) %>%
  summarise(cr.points = mean(cr.points, na.rm = TRUE), .groups = "drop")

xyplot(cr.points ~ prop | CODE, data = ca.avg, type = c("l"),
       group = SPP, auto.key = list(points = FALSE, lines = TRUE, columns = 2))

xyplot(cr.points ~ prop | CODE, 
       data = ca.avg, 
       type = c("l"),
       group = SPP, 
       auto.key = list(points = FALSE, lines = TRUE, columns = 1, title = "Species", space = "right"),
       xlab = "Proportion", 
       ylab = "Crown Points")

write.csv(cr.demog,"Crown_point_summary.csv")


#fit beta distribution for cr.points use ca dataframe (includes each CODE rep df=2)
library(fitdistrplus)
library(agricolae)

beta <- function(ca) {
  ca <- ca[ca$cr.points > 0, ]  
  fit <- fitdist(ca$prop, "beta", method = "mme")  
  return(data.frame(alpha = fit$estimate["shape1"], beta = fit$estimate["shape2"]))
}

beta.fit <- ca %>%
  group_by(BLOCK, CODE, SPP) %>%
  group_modify(~ beta(.x)) %>%
  ungroup()

print(beta.fit)

CODE <- beta.fit %>%
  filter(CODE == "RW", SPP %in% c("RS", "WS"))

#SPP <- beta.fit %>%
  #filter(SPP == "RS" & CODE %in% c("R", "BR", "NR", "RW"))

alpha.aov <- aov(alpha ~ SPP, data = CODE)
beta.aov <- aov(beta ~ SPP, data = CODE)

summary(alpha.aov)
summary(beta.aov)

alpha.tukey <- HSD.test(alpha.aov, "SPP", group = TRUE)
beta.tukey <- HSD.test(beta.aov, "SPP", group = TRUE)

alpha.tukey
beta.tukey

#CODES with significant differences amongst SPP final.ht distributions
#SPP df =1 (2-1); residual df =4 (6-2)

#BN aov alpha p-value: 0.0186 *
#   aov beta p-value: 0.0752 
#   tukey alpha: a b
#   tukey beta: a a

#BR aov alpha p-value: 0.12
#   aov beta p-value: 0.399
#   tukey alpha: a a
#   tukey beta: a a

#BW aov alpha p-value: 0.282
#   aov beta p-value: 0.497
#   tukey alpha: a a
#   tukey beta: a a

#NR aov alpha p-value: 0.525
#   aov beta p-value: 0.335
#   tukey alpha: a a
#   tukey beta: a a

#NW aov alpha p-value: 0.00355 *
#   aov beta p-value: 0.697
#   tukey alpha: a b
#   tukey beta: a a

#RW aov alpha p-value: 0.141
#   aov beta p-value: 0.0184 * 
#   tukey alpha: a a
#   tukey beta: a b


#-------------------------------------------------------------------------------
# Overyielding & Transgressive Overyielding, looking at the plot level
#-------------------------------------------------------------------------------
plot.estimates <- picea %>%
  mutate(qmd.2 = qmd(bapa,tpa),
         rd.2 = bapa/sqrt(qmd.2)) %>%
  group_by(BLOCK, PLOT, CODE) %>%  
  summarise(total.vol = sum(final.vol, na.rm = TRUE),  
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
      
      volume1 <- data$total.vol[data$CODE == species1]
      volume2 <- data$total.vol[data$CODE == species2]
      
      mixture.volume <- data$total.vol[i]
      
      if (length(volume1) > 0 && length(volume2) > 0) {
        average.volume <- mean(c(volume1, volume2), na.rm = TRUE)
        overyielding <- mixture.volume / average.volume
        
        results <- rbind(results, data.frame(Mixture = mixture, Overyielding = overyielding))
      }
    }
  }
  
  return(results)
}

oy.results <- oy(plot.estimates)
print(oy.results)
avg.oy <- oy.results %>%
  group_by(Mixture) %>%
  summarise(avg.oy = mean(Overyielding, na.rm = TRUE))
print(avg.oy)

ggplot(avg.oy, aes(x = Mixture, y = avg.oy)) +
  geom_bar(stat = "identity", fill = "grey") +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  
  labs(x = "Species Mixture", 
       y = NULL) +
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

toy_results <- calculate_transgressive_overyielding(plot.estimates)
avg_toy <- toy_results %>%
  group_by(Mixture) %>%
  summarise(avg_transgressive_oy = mean(Transgressive_Overyielding, na.rm = TRUE))
print(avg_toy)

ggplot(avg_toy, aes(x = Mixture, y = avg_transgressive_oy)) +
  geom_bar(stat = "identity", fill = "grey") +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  
  labs(x = "Species Mixture", 
       y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = c(1))  

#-------------------------------------------------------------------------------
# Multiple Comparisons Test of Stand Level Metrics
# Post-hoc tukey HSD comparison for plot-level metrics
#-------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(agricolae)
library(tibble)

mcp <- picea %>% filter(CODE != "C")

aov <- aov(rdi ~ CODE, data = mcp)

tukey <- HSD.test(aov, "CODE", group = TRUE)

tukey.results <- tukey$groups %>%
  as.data.frame() %>%
  rownames_to_column("CODE") %>%  
  rename(mean_rdi = rdi)

print(tukey.results)

max_rdi <- max(mcp$rdi, na.rm = TRUE)

p.vol.code <- ggplot(mcp, aes(x = factor(CODE), y = rdi)) +
  geom_boxplot(fill = "grey80", color = "black") +  
  geom_text(data = tukey.results, aes(x = CODE, y = max_rdi * 1.05, label = groups), 
            size = 5, fontface = "bold") +  
  ylab("Relative Density Index (%)") +
  xlab("Species Mixture") +
  theme_minimal(base_size = 14)

# Print the plot
print(p.vol.code)

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
    mean.vol = ifelse(CODE == "BR", 164.9, mean.vol),  
    sd.vol = ifelse(CODE == "BR", 137.3, sd.vol)       
  )

ggplot(picea.vol, aes(x = CODE, y = mean.vol, fill = CODE)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean.vol - sd.vol, ymax = mean.vol + sd.vol), width = 0.2) + 
  labs(x = "Species Mixture", y = "Average Plot Volume (ft3/ac)") +
  theme_minimal() +
  theme(legend.position = "none")

#plot.vol with post-hoc tukey HSD comparison results
picea.vol <- left_join(picea.vol, tukey_results, by = "CODE")

ggplot(picea.vol, aes(x = CODE, y = mean.vol, fill = CODE)) +
  geom_bar(stat = "identity", color = "black") +  # Bar plot with border
  geom_errorbar(aes(ymin = mean.vol - sd.vol, ymax = mean.vol + sd.vol), width = 0.2) +  # Error bars
  geom_text(aes(label = groups, y = mean.vol + sd.vol * 1.8), size = 5, fontface = "bold") +  # Tukey labels above bars
  labs(x = "Species Mixture", y = "Plot Volume (ft3/ac)") +
  theme_minimal() +
  theme(legend.position = "none")

# convert to m3/hectare
#picea.vol <- picea %>%
  #filter(CODE != "C") %>% 
  #group_by(CODE) %>%
  #summarise(
    #mean.vol = mean(plot.vol, na.rm = TRUE),
    #sd.vol = sd(plot.vol, na.rm = TRUE)
  #)

#picea.vol <- picea.vol %>%
  #mutate(
    #mean.vol = ifelse(CODE == "BR", 164.9, mean.vol),  
    #sd.vol = ifelse(CODE == "BR", 137.3, sd.vol),
    #mean.vol_m3ha = mean.vol * 0.070,  
    #sd.vol_m3ha = sd.vol * 0.070       
  #)

#ggplot(picea.vol, aes(x = CODE, y = mean.vol_m3ha, fill = CODE)) +
  #geom_bar(stat = "identity") +
  #geom_errorbar(aes(ymin = mean.vol_m3ha - sd.vol_m3ha, ymax = mean.vol_m3ha + sd.vol_m3ha), width = 0.2) + 
  #labs(x = "Species Mixture", y = "Volume (m³/ha)") +
  #theme_minimal() +
  #theme(legend.position = "none")

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

#convert figure to metric 
#picea_converted <- picea %>%
  #mutate(
    #final.ht = final.ht * 0.3048, 
    #DBH.23 = DBH.23 * 2.54         
  #)

# BS plot
#plot_bs <- ggplot(picea_converted %>% filter(SPP == "BS", CODE %in% c("B", "BR", "BW", "BN")) %>%
                    #group_by(CODE, DBH.23) %>%
                    #summarise(avg_final.ht = mean(final.ht, na.rm = TRUE)), 
                  #aes(x = DBH.23, y = avg_final.ht, color = CODE)) +
  #geom_smooth(method = "loess", aes(group = CODE), se = FALSE, linetype = "solid") +
  #labs(x = "DBH (cm)", y = "Height (m)", color = "Species Mixture") +  
  #ggtitle("BS") +
  #theme_minimal() +
  #theme(legend.position = "bottom")

# NS plot
#plot_ns <- ggplot(picea_converted %>% filter(SPP == "NS", CODE %in% c("N", "BN", "NW", "NR")) %>%
                    #group_by(CODE, DBH.23) %>%
                    #summarise(avg_final.ht = mean(final.ht, na.rm = TRUE)), 
                  #aes(x = DBH.23, y = avg_final.ht, color = CODE)) +
  #geom_smooth(method = "loess", aes(group = CODE), se = FALSE, linetype = "solid") +
  #labs(x = "DBH (cm)", y = "Height (m)", color = "Species Mixture") +  
  #ggtitle("NS") +
  #theme_minimal() +
  #theme(legend.position = "bottom")

# RS plot
#plot_rs <- ggplot(picea_converted %>% filter(SPP == "RS", CODE %in% c("R", "BR", "RW", "NR")) %>%
                    #group_by(CODE, DBH.23) %>%
                    #summarise(avg_final.ht = mean(final.ht, na.rm = TRUE)), 
                  #aes(x = DBH.23, y = avg_final.ht, color = CODE)) +
  #geom_smooth(method = "loess", aes(group = CODE), se = FALSE, linetype = "solid") +
  #labs(x = "DBH (cm)", y = "Height (m)", color = "Species Mixture") +  
  #ggtitle("RS") +
  #theme_minimal() +
  #theme(legend.position = "bottom")


# WS plot
#plot_ws <- ggplot(picea_converted %>% filter(SPP == "WS", CODE %in% c("W", "BW", "NW", "RW")) %>%
                    #group_by(CODE, DBH.23) %>%
                    #summarise(avg_final.ht = mean(final.ht, na.rm = TRUE)), 
                  #aes(x = DBH.23, y = avg_final.ht, color = CODE)) +
  #geom_smooth(method = "loess", aes(group = CODE), se = FALSE, linetype = "solid") +
  #labs(x = "DBH (cm)", y = "Height (m)", color = "Species Mixture") +  
  #ggtitle("WS") +
  #theme_minimal() +
  #theme(legend.position = "bottom")

# Combine the plots
#plot_bs + plot_rs + plot_ns + plot_ws + plot_layout(ncol = 2)


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

#convert to metric 
#picea.3 <- picea %>%
  #filter(SPP %in% c("NS", "RS", "WS", "BS"), 
         #CODE != "C",
         #!(CODE == "NR" & SPP == "BS")) %>%  
  #mutate(
    #final.ht = final.ht * 0.3048,  
    #DBH.23 = DBH.23 * 2.54         
  #) %>%
  #group_by(SPP, CODE, DBH.23) %>%
  #summarise(avg.final.ht = mean(final.ht, na.rm = TRUE), .groups = "drop")

#plot.2 <- ggplot(picea.3, aes(x = DBH.23, y = avg.final.ht, color = SPP)) +
  #geom_smooth(method = "loess", aes(group = SPP), se = FALSE, linetype = "solid") +
  #labs(x = "DBH (cm)", y = "Height (m)", color = "Species") +  # Update axis labels
  #theme_minimal() +
  #theme(legend.position = "bottom") +
  #facet_wrap(~CODE, ncol = 2)

#plot.2

#-------------------------------------------------------------------------------
# proportion bapa by the proportion of species in each mixture 
#-------------------------------------------------------------------------------
picea.4 <- picea %>%
  filter(SPP %in% c("NS", "RS", "WS", "BS"), 
         !(CODE %in% c("C", "B", "R", "N", "W")),  
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

#convert to metric
#picea.4 <- picea %>%
  #filter(SPP %in% c("NS", "RS", "WS", "BS"), 
         #!(CODE %in% c("C", "B", "R", "N", "W")),  
         #!(CODE == "NR" & SPP == "BS")) %>%  
  #group_by(SPP, CODE, Proportion) %>%
  #summarise(avg.bapa = mean(bapa, na.rm = TRUE), .groups = "drop") %>%
  #mutate(avg.bapa = avg.bapa * 2.47105)  

#plot.2 <- ggplot(picea.4, aes(x = Proportion, y = avg.bapa, color = SPP)) +
  #geom_smooth(method = "loess", aes(group = SPP), se = FALSE, linetype = "solid") +
  #labs(x = "Proportion", y = "BAPH", color = "Species") +  
  #theme_minimal() +
  #theme(legend.position = "bottom") +
  #facet_wrap(~CODE, ncol = 2)  

#plot.2

#picea.4 <- picea %>%
  #filter(SPP %in% c("NS", "RS", "WS", "BS"), 
         #!(CODE %in% c("C", "B", "R", "N", "W")),  
         #!(CODE == "NR" & SPP == "BS")) %>%  
  #group_by(SPP, CODE, Proportion) %>%
  #summarise(avg.bapa = mean(bapa, na.rm = TRUE), .groups = "drop") %>%
  #mutate(avg.bapa = avg.bapa * 2.47105)  

#plot.2 <- ggplot(picea.4, aes(x = Proportion, y = avg.bapa, color = SPP)) +
  #geom_point(alpha = 0.7, size = 2) +
  #labs(x = "Proportion", y = "BAPH", color = "Species") + 
  #theme_minimal() +
  #theme(legend.position = "bottom") +
  #facet_wrap(~CODE, ncol = 2)  

#plot.2


#-------------------------------------------------------------------------------
# DBH distribution based on monculture vs. mixed Pretzsch and Biber 2016
#-------------------------------------------------------------------------------
spruce.dbh <- picea %>%
  mutate(stand_type = case_when(
    CODE %in% c("NW", "NR", "BN", "BR", "RW", "BW") ~ "Mixed",
    CODE %in% c("N", "W", "B", "R") ~ "Monoculture",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(DBH.23) & !is.na(stand_type))  

ggplot(spruce.dbh, aes(x = DBH.23, fill = stand_type, color = stand_type)) +
  geom_density(alpha = 0.5) +  
  labs(x = "DBH (inches)",
       y = "Density") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(legend.title = element_blank())

ggplot(spruce.dbh, aes(x = DBH.23, fill = stand_type, color = stand_type)) +
  geom_density(alpha = 0.5) +  
  labs(x = "DBH (inches)", y = "Density") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~ CODE)  


spruce.dbh2 <- spruce.dbh %>%
  filter(!(CODE == "BW" & SPP == "RS"),   # Remove SPP RS in CODE BW
         !(CODE == "NR" & SPP == "BS"),   # Remove SPP BS in CODE NR
         !(CODE == "RW" & SPP %in% c("NS", "BS")))  # Remove SPP NS & BS in CODE RW

ggplot(spruce.dbh2, aes(x = DBH.23, fill = SPP, color = SPP)) +
  geom_density(alpha = 0.5) +  
  labs(x = "DBH (in)", y = "Density") +
  scale_x_continuous(breaks = seq(0, max(spruce2$DBH.23, na.rm = TRUE), by = 2),
                     expand = expansion(mult = c(0.05, 0.1))) +  
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~ CODE)

spruce.filtered <- spruce.dbh2 %>% 
  filter(CODE %in% c("BR", "RW", "NR", "R"))

ggplot(spruce.filtered, aes(x = DBH.23, fill = SPP, color = SPP)) +
  geom_density(alpha = 0.5) +  
  labs(x = "DBH (in)", y = "Density", fill = "Species", color = "Species") + 
  scale_x_continuous(breaks = seq(0, max(spruce.filtered$DBH.23, na.rm = TRUE), by = 2),
                     expand = expansion(mult = c(0.05, 0.1))) +  
  theme_minimal() +
  theme(legend.position = "right",  
        legend.title = element_text(size = 10)) +  
  facet_wrap(~ CODE)


#convert to metric 
#ggplot(spruce2, aes(x = DBH.23 * 2.54, fill = SPP, color = SPP)) +
  #geom_density(alpha = 0.5) +  
  #labs(x = "DBH (cm)", y = "Density") +  # Update axis label to cm
  #theme_minimal() +
  #theme(legend.position = "top",
        #legend.title = element_blank()) +
  #facet_wrap(~ CODE)

#-------------------------------------------------------------------------------
# tree-level final.vol distribution based on monculture vs. mixed (Pretzsch and Biber 2016)
#-------------------------------------------------------------------------------
spruce <- picea %>%
  mutate(
    stand_type = case_when(
      CODE %in% c("NW", "NR", "NB", "BR", "RW", "BW", "BN") ~ "Mixed",
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
  geom_density(alpha = 0.5) +  
  labs(x = "Volume (cubic feet)",
       y = "Density") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~ CODE)  

spruce2 <- spruce %>%
  filter(!(CODE == "BW" & SPP == "RS"),   
         !(CODE == "NR" & SPP == "BS"),  
         !(CODE == "RW" & SPP %in% c("NS", "BS")))  

ggplot(spruce2, aes(x = final.vol, fill = SPP, color = SPP)) +
  geom_density(alpha = 0.5) +  
  labs(x = "Volume (ft³)", y = "Density") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~ CODE)

#convert to metric 
#ggplot(spruce2, aes(x = final.vol * 0.0283168, fill = SPP, color = SPP)) +
  #geom_density(alpha = 0.5) +  
  #labs(x = "Volume (m³)", y = "Density") +  # Update axis label to m³
  #theme_minimal() +
  #theme(legend.position = "top",
        #legend.title = element_blank()) +
  #facet_wrap(~ CODE)

#-------------------------------------------------------------------------------
# plot.vol distribution based on monculture vs. mixed (Pretzsch and Biber 2016)
#-------------------------------------------------------------------------------
spruce <- picea %>%
  mutate(
    stand_type = case_when(
      CODE %in% c("NW", "NR", "NB", "BN", "RW", "BW", "BN") ~ "Mixed",
      CODE %in% c("N", "W", "B", "R") ~ "Monoculture",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(plot.vol) & !is.na(stand_type))  

ggplot(spruce, aes(x = plot.vol, fill = stand_type, color = stand_type)) +
  geom_density(alpha = 0.5) +  
  labs(x = "Plot Volume (cubic feet)",
       y = "Density") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(legend.title = element_blank())

ggplot(spruce, aes(x = plot.vol, fill = stand_type, color = stand_type)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  labs(x = "Plot Volume (cubic feet)",
       y = "Density") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~ CODE) 


spruce2 <- spruce %>%
  filter(!(CODE == "BW" & SPP == "RS"),   
         !(CODE == "NR" & SPP == "BS"),  
         !(CODE == "RW" & SPP %in% c("NS", "BS")))  

ggplot(spruce2, aes(x = plot.vol, fill = SPP, color = SPP)) +
  geom_density(alpha = 0.5) +  
  labs(x = "Volume (ft³/ac)", y = "Density") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~ CODE)

#convert to metric 
#ggplot(spruce2, aes(x = plot.vol * 0.070, fill = SPP, color = SPP)) +
  #geom_density(alpha = 0.5) +  
  #labs(x = "Volume (m³/ha)", y = "Density") +  # Update axis label to metric
  #theme_minimal() +
  #theme(legend.position = "top",
        #legend.title = element_blank()) +
  #facet_wrap(~ CODE)
#-------------------------------------------------------------------------------
# final.ht distribution based on monculture vs. mixed (Pretzsch and Biber 2016)
#-------------------------------------------------------------------------------
spruce.ht <- picea %>%
  mutate(
    stand_type = case_when(
      CODE %in% c("NW", "NR", "BN", "BR", "RW", "BW") ~ "Mixed",
      CODE %in% c("N", "W", "B", "R") ~ "Monoculture",
      TRUE ~ NA_character_
    ),
    final.ht = final.ht     
  ) %>%
  filter(!is.na(final.ht) & !is.na(stand_type))  

ggplot(spruce.ht, aes(x = final.ht, fill = stand_type, color = stand_type)) +
  geom_density(alpha = 0.5) +  
  labs(x = "Height (m)",  
       y = "Density") +
  scale_fill_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  scale_color_manual(values = c("Mixed" = "blue", "Monoculture" = "red")) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~ CODE)

spruce.ht2 <- spruce.ht %>%
  filter(!(CODE == "BW" & SPP == "RS"),   
         !(CODE == "NR" & SPP == "BS"),  
         !(CODE == "RW" & SPP %in% c("NS", "BS")))  

ggplot(spruce.ht2, aes(x = final.ht, fill = SPP, color = SPP)) +
  geom_density(alpha = 0.5) +  
  labs(x = "Height (ft)", y = "Density") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~ CODE) #weibull looks good to fit these, if skewed consider a gamma or log-normal distribution

#convert to metric
#ggplot(spruce2, aes(x = final.ht * 0.3048, fill = SPP, color = SPP)) +
  #geom_density(alpha = 0.5) +  
  #labs(x = "Height (m)", y = "Density") +  # Update axis label to meters
  #theme_minimal() +
  #theme(legend.position = "top",
        #legend.title = element_blank()) +
  #facet_wrap(~ CODE)

#-------------------------------------------------------------------------------
# fitting Weibull distribution to DBH.23 distrubution
#-------------------------------------------------------------------------------
# account for replicates 
# [1] calculate the weibull shape/scale for each block/plot/code/species, and just use a mcp for both shape and scale by group

library(fitdistrplus) 
library(multcomp) 

spruce.filt <- picea %>%
  filter(!((CODE == "BW" & SPP == "RS") | 
             (CODE == "NR" & SPP == "BS") | 
             (CODE == "RW" & SPP %in% c("NS", "BS"))) & 
           !is.na(DBH.23) & !is.nan(DBH.23))

fitted.mods2 <- spruce.filt %>%
  split(~ BLOCK + CODE + SPP) %>%
  map(~ tryCatch(fitdist(.x$DBH.23, "weibull"), error = function(e) NULL)) %>%
  keep(~ !is.null(.)) 
fitted.mods2

# extract shape and scale parameters into a tidy dataframe
ss.df <- fitted.mods2 %>%
  imap_dfr(~ tibble(
    BLOCK = strsplit(.y, "_")[[1]][1],
    #PLOT = strsplit(.y, "_")[[1]][2],
    CODE = strsplit(.y, "_")[[1]][3],
    SPP = strsplit(.y, "_")[[1]][4],
    Shape = .x$estimate[["shape"]],
    Scale = .x$estimate[["scale"]]
  ))


ss.df <- ss.df %>%
  mutate(BLOCK = as.character(BLOCK)) %>%
  separate(BLOCK, into = c("BLOCK", "CODE", "SPP"), sep = "\\.", remove = FALSE)

#ss.df.avg <- ss.df %>%
  #group_by(BLOCK + CODE, SPP) %>%
  #summarize(
    #Shape = mean(Shape, na.rm = TRUE),
    #Scale = mean(Scale, na.rm = TRUE),
    #.groups = "drop"  
  #)

#head(ss.df.avg)

#tukey <- function(ss.df) {
  #mod.shape <- aov(Shape ~ CODE + SPP + BLOCK, data = ss.df)
  #mod.scale <- aov(Scale ~ CODE + SPP + BLOCK, data = ss.df)
  
  #tukey.shape <- TukeyHSD(mod.shape)
  #tukey.scale <- TukeyHSD(mod.scale)
  
  #list(Tukey.Shape = tukey.shape, Tukey.Scale = tukey.scale)
#}

#t.results <- ss.df %>%
  #group_by(CODE, SPP) %>%  
  #group_split() %>%        
  #map(tukey)

#t.results


# test differences between SPP Scale and Shape of DBH distributions 1 code at a time
CODE <- ss.df %>%
  filter(CODE == "BR", SPP %in% c("BS", "RS"))

SPP <- ss.df %>%
  filter(SPP == "WS" & CODE %in% c("W", "BW", "NW", "RW"))

# ANOVA and Tukey HSD test
shape.aov <- aov(Shape ~ CODE, data = SPP)
scale.aov <- aov(Scale ~ CODE, data = SPP)

summary(shape.aov)
summary(scale.aov)

shape.tukey <- HSD.test(shape.aov, "CODE", group = TRUE)
scale.tukey <- HSD.test(scale.aov, "CODE", group = TRUE)

shape.tukey
scale.tukey


#CODES with significant differences amongst SPP DBH.23 distributions
#SPP df =1 (2-1); residual df =4 (6-2)

#BN aov shape p-value: 0.392
#   aov scale p-value: 0.112 
#   tukey shape: a a 
#   tukey scale: a a
#BR aov shape p-value: 0.488
#   aov scale p-value: 0.479
#   tukey shape: a a 
#   tukey scale: 
#BW aov shape p-value: 0.237
#   aov scale p-value: 0.734
#   tukey shape: a a
#   tukey scale: a a
#NR aov shape p-value: 0.00334
#   aov scale p-value: 0.000535
#   tukey shape: a b
#   tukey scale: a b
#NW aov shape p-value: 0.0815
#   aov scale p-value: 0.367
#   tukey shape: a a
#   tukey scale: a a
#RW aov shape p-value: 0.0155
#   aov scale p-value: 4e-04
#   tukey shape: a b
#   tukey scale: a b

#SPP with significant differences amongst DBH.23 distributions in CODE
# CODE df = 4; residual df = 8

#BS in B, BN, BR, BW
# no significant differences at 0.05 in shape & scale

#NS in N, BN, NR, NW
# no significant differences at 0.05 in shape & scale

#RS in R, BR, NR, RW
# no significant differences at 0.05 in shape & scale

#WS in W, BW, NW, RW
# no significant differences at 0.05 in shape & scale



# student t-test
#shape.ttest <- t.test(Shape ~ SPP, data = CODE)

#scale.ttest <- t.test(Scale ~ SPP, data = CODE)

#shape.ttest
#scale.ttest




# [2] use 1-inch dbh classes and use this your response variable, fit a weibull or other distribution with nlme approach. 

dbh1 <- spruce2 %>%
  mutate(DBH1 = cut(DBH.23, breaks = seq(0, max(DBH.23), by = 1), right = FALSE, labels = FALSE)) %>%
  group_by(SPP, CODE) %>%
  ungroup()

dbh1 <- dbh1[!is.na(dbh1$DBH1), ]

weibull_model <- function(x, shape, scale) {
  return(scale * (shape - 1) * (x / scale)^(shape - 1) * exp(-(x / scale)^shape))
}

weibull.mod <-  nlme(DBH1 ~  weibull_model(DBH1, shape, scale),
                data = dbh1,
                fixed =shape + scale ~ CODE * SPP,
                random = shape + scale ~ CODE * SPP| BLOCK,  
                na.action = na.pass,
                start = c(shape = 1, scale = 1),
                control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))


#mixed effects weibull model
#spruce2 <- spruce2 %>%
  #filter(!is.na(DBH.23) & DBH.23 > 0)  
#spruce.2 <- spruce2 %>% mutate(log.DBH = log(DBH.23))  # log-transform for Weibull-like fit

#model.mixed <- lmer(log.DBH ~ CODE + (1 | BLOCK), data = spruce.2)
#summary(model.mixed)
#model.null <- lmer(log_DBH ~ (1 | BLOCK), data = spruce.2)  #
#anova(model.null, model.mixed)


#-------------------------------------------------------------------------------
# fitting Weibull distribution to final.ht distributions
#-------------------------------------------------------------------------------
# account for replicates 
# [1] calculate the weibull shape/scale for each block/plot/code/species, and just use a mcp for both shape and scale by group

library(fitdistrplus) 
library(multcomp) 

height.filt <- picea %>%
  filter(!((CODE == "BW" & SPP == "RS") | 
             (CODE == "NR" & SPP == "BS") | 
             (CODE == "RW" & SPP %in% c("NS", "BS"))) & 
           !is.na(final.ht) & !is.nan(final.ht))

fitted.mods3 <- height.filt %>%
  split(~ BLOCK + CODE + SPP) %>%
  map(~ tryCatch(fitdist(.x$final.ht, "weibull"), error = function(e) NULL)) %>%
  keep(~ !is.null(.)) 
fitted.mods3

# extract shape and scale parameters into a tidy dataframe
ss.ht.df <- fitted.mods3 %>%
  imap_dfr(~ tibble(
    BLOCK = strsplit(.y, "_")[[1]][1],
    #PLOT = strsplit(.y, "_")[[1]][2],
    CODE = strsplit(.y, "_")[[1]][3],
    SPP = strsplit(.y, "_")[[1]][4],
    Shape = .x$estimate[["shape"]],
    Scale = .x$estimate[["scale"]]
  ))


ss.ht.df <- ss.ht.df %>%
  mutate(BLOCK = as.character(BLOCK)) %>%
  separate(BLOCK, into = c("BLOCK", "CODE", "SPP"), sep = "\\.", remove = FALSE)

# test differences between SPP Scale and Shape of final.ht distributions 1 code at a time
CODE <- ss.ht.df %>%
  filter(CODE == "RW", SPP %in% c("RS", "WS"))

#SPP <- ss.ht.df %>%
  #filter(SPP == "BS" & CODE %in% c("B", "BN", "BR", "BW"))


# ANOVA and Tukey HSD test
shape.aov <- aov(Shape ~ SPP, data = CODE)
scale.aov <- aov(Scale ~ SPP, data = CODE)

summary(shape.aov)
summary(scale.aov)

shape.tukey <- HSD.test(shape.aov, "SPP", group = TRUE)
scale.tukey <- HSD.test(scale.aov, "SPP", group = TRUE)

shape.tukey
scale.tukey


#CODES with significant differences amongst SPP final.ht distributions
#SPP df =1 (2-1); residual df =4 (6-2)

#BN aov shape p-value: 0.971
#   aov scale p-value: 0.149 
#   tukey shape: a a 
#   tukey scale: a a
#BR aov shape p-value: 0.333
#   aov scale p-value: 0.0796
#   tukey shape: a a 
#   tukey scale: a a 
#BW aov shape p-value: 0.0821
#   aov scale p-value: 0.6
#   tukey shape: a a
#   tukey scale: a a
#NR aov shape p-value: 0.0843
#   aov scale p-value: 0.238
#   tukey shape: a a
#   tukey scale: a a
#NW aov shape p-value: 0.991
#   aov scale p-value: 0.172
#   tukey shape: a a
#   tukey scale: a a
#RW aov shape p-value: 0.112
#   aov scale p-value: 0.00114 *
#   tukey shape: a a
#   tukey scale: a b

#SPP with significant differences amongst final.ht distributions in CODE
# CODE df = 4; residual df = 8

#BS in B, BN, BR, BW
#   shape p-value 0.719
#   B a, BN a, BR a, BW a
#   scale p-value: 0.864
#   B a, BN a, BR a, BW a

#NS in N, BN, NR, NW
#   shape p-value 0.535
#   N a, BN a, NR a, NW a
#   scale p-value: 0.57
#   N a, BN a, NR a, NW a

#RS in R, BR, NR, RW
#   shape p-value 0.864
#   R a, BR a, NR a, RW a
#   scale p-value: 0.326
#   R a, BR a, NR a, RW a

#WS in W, BW, NW, RW
#   shape p-value 0.0319 *
#   RW a, NW ab, BW b, W b
#   scale p-value: 0.74
#   RW a, NW a, BW a, W a