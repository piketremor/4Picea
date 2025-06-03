#install.packages("devtools")
#install.packages("G:/My Drive/MEForLab",
#                 repos = NULL,
#                 type = "source")
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
library(tidyr)


#setwd("C:/Home/Documents/GitHub/4Picea/code/raw")

setwd("G:/Shared drives/4Picea/4Picea/raw")
setwd("~/Documents/GitHub/4Picea/code/raw")

picea <- read.csv("4Picea.csv")
stemform <- read.csv("StemForm.csv")
stemform <- stemform %>%
  mutate(TREE = as.character(TREE)) %>%
  distinct(BLOCK, PLOT, TREE, .keep_all = TRUE) #duplicates in stem form data set by TREE
site <- read.csv("4Picea_30m.csv")

head(picea)
names(picea)
#-------------------------------------------------------------------------------
#join stemform to picea by BLOCK, PLOT, TREE #join site to picea by BLOCK, PLOT
#-------------------------------------------------------------------------------
picea <- picea %>%
  mutate(TREE = as.character(TREE)) %>%
  left_join(stemform %>% mutate(TREE = as.character(TREE)), by = c("BLOCK", "PLOT", "TREE", "SPP")) %>%
  left_join(site, by = c("BLOCK", "PLOT")) %>%
  filter(SPP %in% c("BS", "NS", "RS", "WS"), CODE != "C")

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
species.prop %>% print(n = Inf)

#library(writexl)
#write_xlsx(species.prop, path = "species_prop.xlsx")

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

#bapa.tpa.summary %>% print(n = Inf)
#write_xlsx(bapa.tpa.summary, path = "bapa.tpa.plot.summary.xlsx")

#BR.dbh <- picea %>%
#  filter(CODE == "BR") %>%
#  group_by(BLOCK, PLOT) %>%
#  summarise(
#    min_DBH = min(DBH.23, na.rm = TRUE),
#    max_DBH = max(DBH.23, na.rm = TRUE),
#    .groups = 'drop'
#  )

#print(BR.dbh)

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

#-------------------------------------------------------------------------------
# #DBH distribution
#-------------------------------------------------------------------------------
#picea$dbh.class <- 2*as.integer((picea$DBH.23+(2/2))/2)

#picea.dist <- picea %>%
#  filter(STATUS.23 == "1") %>%
#  mutate(DiameterClass = cut(DBH.23, 
#                             breaks = seq(0, 15, by = 1), 
#                             right = FALSE)) %>%
#  group_by(dbh.class) %>%
#  summarise(Count = n(), 
#            .groups = 'drop')

#ggplot(picea.dist, aes(x = dbh.class, y = Count, fill = dbh.class)) +
#  geom_bar(stat = "identity", show.legend = FALSE) +
#  labs(title = "Diameter Class Distribution",
#       x = "Diameter Classes (DBH in inches)",
#       y = "Number of Trees") +
#  scale_x_discrete(drop = FALSE) +
#  theme_minimal() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)

# Diameter distribution by CODE
#picea$dbh.class <- 2 * as.integer((picea$DBH.23 + (2 / 2)) / 2)

#picea.dist <- picea %>%
#  filter(STATUS.23 == "1") %>%
#  mutate(DiameterClass = cut(DBH.23, 
#                             breaks = seq(0, max(picea$DBH.23, na.rm = TRUE), by = 2), 
#                             right = FALSE)) %>%
#  group_by(CODE, DiameterClass) %>%
#  summarise(Count = n(), 
#            .groups = 'drop')

#ggplot(picea.dist, aes(x = DiameterClass, y = Count, fill = DiameterClass)) +
#  geom_bar(stat = "identity", show.legend = FALSE) +
#  facet_wrap(~ CODE) +  # Create separate plots for each CODE
#  labs(x = "Diameter Class (DBH in inches)",
#       y = "Number of Trees") +
#  scale_x_discrete(labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12"), drop = FALSE) +  
#  theme_minimal() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for better readability

#-------------------------------------------------------------------------------
# Generating height model
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
                random = a + b ~ 1 | SPP,  
                na.action = na.pass,
                start = c(a = 4.5, b = -6),
                control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))
performance(ht.mod)
summary(ht.mod)
ranef(ht.mod)
AIC(ht.mod) #7029

predictions <- fitted(ht.mod)
observed <- spruce$HT.23

observed <- observed[!is.na(observed) & !is.na(predictions)]
predictions <- predictions[!is.na(observed) & !is.na(predictions)]

mb <- mean(observed - predictions)
mab <- mean(abs(observed - predictions))

sst <- sum((observed - mean(observed))^2) 
sse <- sum((observed - predictions)^2)  
rsq <- 1 - sse / sst  

n <- length(observed)  
p <- length(fixed.effects(ht.mod)) 
adjrsq <- 1 - ((1 - rsq) * (n - 1) / (n - p - 1))  

mb
mab
rsq
adjrsq


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


#-------------------------------------------------------------------------------
# HDR - (tree level)
#-------------------------------------------------------------------------------
hdr <- picea %>%
  mutate(
    DBH.23.m = DBH.23 * 0.0254, #convert DBH.23 from inches to meters
    final.ht.m = final.ht * 0.3048, #convert final.ht from ft to meters
    HDR = final.ht.m / DBH.23.m
  ) %>%
  filter(is.finite(HDR)) %>%
  group_by(CODE, BLOCK, PLOT, SPP) %>%
  summarise(
    mean_HDR_replicate = mean(HDR, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  group_by(CODE, SPP) %>%
  summarise(
    mean_HDR = mean(mean_HDR_replicate, na.rm = TRUE),
    sd_HDR   = sd(mean_HDR_replicate, na.rm = TRUE),
    min_HDR  = min(mean_HDR_replicate, na.rm = TRUE),
    max_HDR  = max(mean_HDR_replicate, na.rm = TRUE),
    .groups = 'drop'
  )

print(hdr, n = Inf)

#-------------------------------------------------------------------------------
# Basal area larger (bal) - (tree level)
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
       subset = CODE %in% c("R", "RW", "W", "N", "NR", "NW", "B", "BN", "BR", "BW"),
       layout = c(2, ceiling(length(unique(spruce.bal$CODE))/2))) 

picea <- picea %>%
  left_join(spruce.bal %>% select(uid, BLOCK, PLOT, TREE, bal), 
            by = c("uid", "BLOCK", "PLOT", "TREE"))

#-------------------------------------------------------------------------------
#Height larger (htl) - (tree level)
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
# Maximum crown width (tree level)
#-------------------------------------------------------------------------------
picea <- picea %>%
  mutate(DBH = as.numeric(DBH.23), 
         MCW = mapply(MCW, SPP = SPP, DBH = DBH))

#-------------------------------------------------------------------------------
# CCF 
#-------------------------------------------------------------------------------
spruce.ccf <- picea%>%
  mutate(CW = mapply(MCW,SPP=SPP,DBH=DBH),
         CA = (CW/2)^2*pi,
         CA.exp = CA*10)%>%
  group_by(BLOCK,PLOT)%>%
  summarize(CCF = (sum(CA.exp)/43560)*100)
  
picea <- picea %>%
  left_join(spruce.ccf %>% select(BLOCK, PLOT, CCF), 
            by = c("BLOCK", "PLOT"))

--------------------------------------------------------------------
# Site Index - (calculated at tree level, summarize by Block and Plot)
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
# Top Height/HT40 (vicary height) - (tree or plot level?)
#-------------------------------------------------------------------------------
picea <- picea%>%
  mutate(topht=mapply(vicary.height,SPP="RS", age=28, si=vicary.si))

picea %>%
  group_by(BLOCK, PLOT) %>%
  summarise(spruce.topht = mean(topht, na.rm = TRUE))

#-------------------------------------------------------------------------------
# qmd - (plot level)
#-------------------------------------------------------------------------------
spruce.qmd <- picea %>%
  group_by(BLOCK, PLOT) %>%
  summarise(qmd = sqrt(sum(DBH.23^2, na.rm = TRUE) / n()), .groups = "drop")

mean(spruce.qmd$qmd)

picea <- picea %>%
  left_join(spruce.qmd %>% select(BLOCK, PLOT, qmd), 
            by = c("BLOCK", "PLOT"))

#-------------------------------------------------------------------------------
# Relative density index - (plot level)
#-------------------------------------------------------------------------------
spruce.rd <- picea %>%
  mutate(rdi = mapply(relative.density.index, bapa, qmd))

spruce.rd.summary <- spruce.rd %>%
  group_by(BLOCK, PLOT) %>%
  summarise(rdi = mean(rdi, na.rm = TRUE))  # Ensure only one value per BLOCK-PLOT

picea <- picea %>%
  left_join(spruce.rd.summary, by = c("BLOCK", "PLOT"))


#-------------------------------------------------------------------------------
# Relative density index by SPP within each PLOT
#-------------------------------------------------------------------------------
picea <- picea %>%
  mutate(rdi = mapply(relative.density.index, bapa, qmd))

# summarize RDI by PLOT and SPP
rd.species.plot <- picea %>%
  group_by(BLOCK, PLOT, SPP) %>%
  summarise(rd.spp = sum(rdi, na.rm = TRUE), .groups = "drop")

# calculate total rdi by PLOT
rd.plot <- rd.species.plot %>%
  group_by(BLOCK, PLOT) %>%
  summarise(rd.total = sum(rd.spp, na.rm = TRUE), .groups = "drop")

#join and calculate RD % by SPP
rd <- rd.species.plot %>%
  left_join(rd.plot, by = c("BLOCK", "PLOT")) %>%
  mutate(rd = rd.spp / rd.total) %>%
  select(BLOCK, PLOT, SPP, rd)


rd.1 <- rd %>% #create SPP columns
  pivot_wider(names_from = SPP, values_from = rd, names_prefix = "", values_fill = 0) %>%
  rename_with(.fn = ~ paste0(.x, ".rd"), .cols = c("BS", "NS", "RS", "WS"))

site <- site %>%
  left_join(rd.1, by = c("BLOCK", "PLOT"))

#-------------------------------------------------------------------------------
# Relative spacing - (plot level)
#-------------------------------------------------------------------------------
str(relative.spacing)

picea <- picea %>%
  mutate(rs=mapply(relative.spacing,tpa=tpa, domht=topht))

picea %>%
  group_by(BLOCK, PLOT) %>%
  summarise_at(vars(rs), list(name = mean))

#-------------------------------------------------------------------------------
# Stand density index - (plot level)
#-------------------------------------------------------------------------------
picea <- picea %>%
  mutate(sdi=mapply(stand.density.index,tpa=tpa, qmd=qmd))

picea %>%
  group_by(BLOCK, PLOT) %>%
  summarise_at(vars(sdi), list(name = mean))

xyplot(qmd~tpa|CODE,data=picea)

#-------------------------------------------------------------------------------
# Honer volume calculation at the individual tree level (ft3)
#-------------------------------------------------------------------------------
#picea <- picea%>%
  #mutate(vol=mapply(vol.calc,SPP=SPP,DBH=DBH.23,HT=final.ht))
#xyplot(vol~DBH.23|CODE,data=picea)

#calculate volume at plot level
#picea <- picea %>%
  #group_by(BLOCK, PLOT) %>%
  #mutate(stand.vol = sum(vol, na.rm = TRUE) * 10)

#plot(vol2$fitted.values, resid(vol), 
     #xlab = "Fitted Values", 
     #ylab = "Residuals", 
     #main = "Residual Plot",
     #pch = 20, col = "black")
#abline(h = 0, lty = 2, col = "red")

#-------------------------------------------------------------------------------
# Kozak volume calculation at the individual tree level (ft3)
#-------------------------------------------------------------------------------
require(devtools)
devtools::install("GreenTimbMerch")
require(GreenTimbMerch)

picea$dbh.cm <- picea$DBH.23*2.54
picea$ht.m <- picea$ht.fit/3.28
spruce <- picea %>%
  dplyr::filter(SPP %in% c("NS", "WS", "RS", "BS"),
                CODE != "C")


spruce$new.vol.m3 <- mapply(KozakTreeVol,'ib',spruce$SPP,spruce$dbh.cm,spruce$ht.m,
                            Planted=TRUE)
spruce <- spruce %>%
  group_by(BLOCK, PLOT) %>%
  mutate(
    p.vol.ha = sum(new.vol.m3, na.rm = TRUE) * 10 * 2.47105
  )

spruce <- spruce %>%
  group_by(BLOCK, PLOT) %>%
  mutate(
    p.vol.ac = p.vol.ha * 14.28
  )

spruce.vol <- spruce %>%
  group_by(BLOCK, PLOT) %>%
  summarise(
    p.vol.ha = sum(new.vol.m3, na.rm = TRUE) * 10 * 2.47105,
    .groups = "drop"
  ) %>%
  mutate(
    p.vol.ac = p.vol.ha * 14.28
  )

picea <- picea %>%
  left_join(spruce.vol, by = c("BLOCK", "PLOT"))


#------------------------------------------------------------------------------
# Relative volume (p.vol.ha) by SPP within each PLOT
#------------------------------------------------------------------------------
#summary of volume per species per plot
vol.spp.plot <- spruce %>%
  group_by(BLOCK, PLOT, SPP) %>%
  summarise(vol.spp = sum(p.vol.ha, na.rm = TRUE), .groups = "drop")

# total plot volume
vol.plot <- vol.spp.plot %>%
  group_by(BLOCK, PLOT) %>%
  summarise(vol.total = sum(vol.spp, na.rm = TRUE), .groups = "drop")

# volume percent by SPP
vol <- vol.spp.plot %>%
  left_join(vol.plot, by = c("BLOCK", "PLOT")) %>%
  mutate(vol.rel = vol.spp / vol.total) %>%
  select(BLOCK, PLOT, SPP, vol.rel)

vol.1 <- vol %>%
  pivot_wider(names_from = SPP, values_from = vol.rel, names_prefix = "", values_fill = 0) %>%
  rename_with(.fn = ~ paste0(.x, ".vol"), .cols = c("BS", "NS", "RS", "WS"))

site <- site %>%
  left_join(vol.1, by = c("BLOCK", "PLOT"))

#-------------------------------------------------------------------------------
# Variable selection procedures using VSURF package for integer response variable
#-------------------------------------------------------------------------------
require(randomForest)
require(VSURF)
require(caTools)
require(pscl)

spruce$deductclass <- 5*as.integer((spruce$T_DEDUCT+(5/2))/5)
hist(spruce$deductclass,
     xlab = "Deduction (%)",
     ylab = "Number of Trees",
     main = "Distribution of Deduction Classes")
d.set <- spruce
d.set$T_DEDUCT[is.na(d.set$T_DEDUCT)] <- 999
d.set <- dplyr::filter(d.set,T_DEDUCT<998)
#d.set$HT.23[is.na(d.set$HT.23)] <- 0
#d.set <- dplyr::filter(d.set,HT.23>0)
plot(d.set$DBH.23,d.set$HT.23)

names(d.set)
#obs <- d.set$deductclass
#obs[is.na(obs)] <- 0
#preds <- d.set[c(4,11,12,21:46,51:54,57:61,65,66,71,76:79)]
#preds[is.na(preds)] <- 0
#vs <- VSURF(preds,obs,ncores = 4)
#vs$varselect.pred
#names(preds)

# SPP, CODE

#preds2 <- d.set[c(11,12,21:46,51:55,57:61,65,66,71,76:79)]
#preds2[is.na(preds2)] <- 0
#obs <- d.set$deductclass
#obs[is.na(obs)] <- 0
#nvs <- VSURF(preds2,obs)
#nvs$varselect.pred
#names(preds2)

#CODE, topht, bapa, bal


mod3.1 <- zeroinfl(deductclass~rs+sdi+bal+roughness+CODE|sdi+bal+CODE+SPP,
                 data=d.set,dist="negbin")
summary(mod3.1)
performance(mod3.1)

#mod3.2 <- zeroinfl(deductclass~topht+bapa+bal+CODE|topht+bapa+bal+CODE+SPP,
                   #data=d.set,dist="negbin")
#summary(mod3.2)
#AIC(mod3.2, mod3.1)

predictions <- predict(mod3.1, type = "response")
observed <- d.set$deductclass
mb <- mean(observed - predictions)
mab <- mean(abs(observed - predictions))
mb
mab

# dang, that is good. 

spruce$fit.deduction <- ifelse(spruce$SPP=="RS"|
                              spruce$SPP=="WS"|
                              spruce$SPP=="BS"|
                              spruce$SPP=="NS",predict(mod3.1,type="response"),0)
spruce$fit.deduct.class <- (5*as.integer((spruce$fit.deduction+(5/2))/5))/100
# will use the modeled reduction
spruce$fit.deduct.class[spruce$fit.deduct.class>1] <- 1

spruce$adj.vol <- spruce$p.vol.ac*spruce$fit.deduct.class
plot(spruce$p.vol.ac,spruce$adj.vol,ylim=c(0,2000),xlim=c(0,2000))
spruce$final.vol <- spruce$p.vol.ac-spruce$adj.vol

#convert final.vol to m3/ha
spruce <- spruce %>%
  mutate(merch.vol.ha = final.vol * 0.0698)


vol.summary <- spruce %>%
  group_by(BLOCK, PLOT,CODE) %>%
  summarise(
    new_vol_m3 = mean(new.vol.m3, na.rm = TRUE),
    p_vol_ha = mean(p.vol.ha, na.rm = TRUE),
    p_vol_ac = mean(p.vol.ac, na.rm = TRUE),
    deduct_class = mean(deductclass, na.rm = TRUE),
    fit_deduction = mean(fit.deduction, na.rm = TRUE),
    fit_deduct_class = mean(fit.deduct.class, na.rm = TRUE),
    adj_vol = mean(adj.vol, na.rm = TRUE),
    final_vol = mean(final.vol, na.rm = TRUE),
    merch_vol_ha = mean(merch.vol.ha, na.rm = TRUE),
    .groups = "drop"
  )

#write_xlsx(vol.summary, path = "volume.plot.summary.xlsx")


# Test for significant differences amongst p.vol.ha and merch.vol.ha by CODE
library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)

pie <- spruce %>%
  dplyr::select(CODE, p.vol.ha, merch.vol.ha) %>%
  pivot_longer(cols = c(p.vol.ha, merch.vol.ha),
               names_to = "VolumeType",
               values_to = "Volume")
pie$VolumeType <- as.factor(pie$VolumeType)

# run two-way ANOVA with interaction
anova.model <- aov(Volume ~ CODE * VolumeType, data = pie)
summary(anova.model)

# post-hoc pairwise comparisons
library(emmeans)
emm <- emmeans(anova.model, ~ VolumeType | CODE)
pairwise.comparisons <- contrast(emm, method = "pairwise", adjust = "tukey")
summary(pairwise.comparisons)

#summarize for graphic
volume.summary <- pie %>%
  group_by(CODE, VolumeType) %>%
  summarise(mean.volume = mean(Volume, na.rm = TRUE),
            se.volume = sd(Volume, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

volume.summary$VolumeType <- factor(volume.summary$VolumeType, levels = c("p.vol.ha", "merch.vol.ha"))

ggplot(volume.summary, aes(x = CODE, y = mean.volume, fill = VolumeType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean.volume - se.volume, ymax = mean.volume + se.volume),
                position = position_dodge(width = 0.9),
                width = 0.2) +
  scale_fill_manual(values = c("p.vol.ha" = "grey", "merch.vol.ha" = "blue"),
                    labels = c("p.vol.ha" = "Volume", "merch.vol.ha" = "Merchantable Volume")) +  # Set custom labels
  labs(x = "Treatment",
       y = "Volume (m³ ha⁻¹)",
       fill = "Volume Type") +
  theme_minimal()


#-------------------------------------------------------------------------------
# HCB Model
#-------------------------------------------------------------------------------
#h.set <- picea
#h.set$HCB.23[is.na(h.set$HCB.23)] <- 999
#h.set <- dplyr::filter(h.set,HCB.23<998)

#names(h.set)
#obs <- h.set$HCB.23
#obs[is.na(obs)] <- 0
#preds <- h.set[c(4,6,7,11,12,21:46,50:61,65,66,71,72,75,76,77,78,79,80,81,82)]
#preds2 <- h.set[c(4,6,7,11,12,21:46,50:54,57:61)] #take out stand structure metrics
#preds[is.na(preds2)] <- 0
#vs <- VSURF(preds2,obs,ncores = 4)
#vs$varselect.pred
#names(preds2)

#HT.23, log(bal), log(CCF), hd, DBH.23 for A. Weiskittel paper
#SPP, CODE, rdi, vicary.si (circular), bal for preds
#SPP, HT.23, DBH.23, CODE for preds2
#SPP, HT.23, DBH.23, vpdmax

hcb.mod1 <- lm(HCB.23 ~final.ht + DBH.23  + vpdmax + factor(SPP), data = picea) 
summary(hcb.mod1)
performance(hcb.mod1)

hcb.mod4 <- lm(HCB.23 ~final.ht + DBH.23 + factor(SPP) + factor(CODE),data = picea) #log final.ht or DBH.23 did not improve AIC
summary(hcb.mod4)
performance(hcb.mod4)
AIC(hcb.mod1, hcb.mod4) #hcb.mod4 is better model, hcb.mod1 had a higher R2 but also high AIC 

rmse <- rmse(hcb.mod4)
print(paste("RMSE:", rmse))

MB <- mean(residuals(hcb.mod4))
print(paste("Mean Bias (MB):", MB))

MAB <- mean(abs(residuals(hcb.mod4)))
print(paste("Mean Absolute Bias (MAB):", MAB))

# predict hcb.mod4 on all spruce, and join to picea df
picea$hcb.fit <- ifelse(picea$SPP=="RS"|
                          picea$SPP=="WS"|
                          picea$SPP=="BS"|
                          picea$SPP=="NS",predict(hcb.mod4,type="response"),0)
picea$final.hcb <- ifelse(!is.na(picea$HCB.23), picea$HCB.23, picea$hcb.fit)

#-------------------------------------------------------------------------------
# Crown Stratification
#-------------------------------------------------------------------------------
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
al <- spruce.only[c(1:8,12,66:80,87)]
head(al)
pa <- left_join(nd2,al)
head(pa)
pa$top.prop <- (pa$final.ht-pa$Low)/(pa$Top-pa$Low)
pa$bot.prop <- (pa$fit.hcb-pa$Low)/(pa$Top-pa$Low)
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

#write.csv(cr.demog,"Crown_point_summary.csv")

ca.avg.smooth <- ca.avg %>%
  group_by(CODE, SPP) %>%
  do({
    smoothed <- smooth.spline(.$prop, .$cr.points, spar = 0.7)
    data.frame(prop = smoothed$x, cr.points = pmax(smoothed$y, 0))  # Ensures no negative values
  }) %>%
  ungroup()

xyplot(prop~cr.points | CODE, 
       data = ca.avg.smooth, 
       type = "l", 
       group = SPP, 
       auto.key = list(points = FALSE, lines = TRUE, columns = 1, title = "Species", space = "right"),
       xlab = "Crown Points", 
       ylab = "Proportion")

# Area under the curve: to determine is stratification influences LAI or vol
library(car)

interval <- ca.avg.smooth %>%
  group_by(CODE, SPP) %>%
  summarise(
    min_prop = round(min(prop, na.rm = TRUE), 1), 
    max_prop = round(max(prop, na.rm = TRUE), 1), 
    range_prop = round(max_prop - min_prop, 1),  
    .groups = "drop"
  )
interval

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

alpha.aov <- aov(alpha ~ SPP, data = CODE)
beta.aov <- aov(beta ~ SPP, data = CODE)

summary(alpha.aov)
summary(beta.aov)

#CODES with significant differences amongst SPP CPI
#SPP df =1 (2-1); residual df =4 (6-2)

#BN aov alpha p-value: 0.00054 ***
#   aov beta p-value: 0.0667

#BR aov alpha p-value: 0.0489*
#   aov beta p-value: 0.48

#BW aov alpha p-value: 0.416
#   aov beta p-value: 0.493

#NR aov alpha p-value: 0.339
#   aov beta p-value: 0.282

#NW aov alpha p-value: 0.00464 **
#   aov beta p-value: 0.173

#RW aov alpha p-value: 0.00412 **
#   aov beta p-value: 0.0184 * 


library(multcomp)

beta.fit$SPP <- as.factor(beta.fit$SPP)
beta.fit$CODE <- as.factor(beta.fit$CODE)

amod <- lm(alpha ~ SPP + CODE, data = beta.fit)  
bmod <- lm(beta ~ SPP + CODE, data = beta.fit)  

summary(amod)
summary(bmod)

alpha.glht <- glht(amod, linfct = mcp(CODE = "Tukey")) 
summary(alpha.glht, p.adjust.method = "bonferroni")

beta.glht <- glht(bmod, linfct = mcp(CODE = "Tukey"))
summary(beta.glht, p.adjust.method = "bonferroni") 

#-------------------------------------------------------------------------------
# Volume predictions for all treatments
#-------------------------------------------------------------------------------
plot.summary <- spruce %>%
  mutate(
    ba = DBH.23^2 * 0.005454,  
    ef = 10,
    c.area = (MCW/2)^2*pi,
    ht.m = ht.fit/3.28,
    dbh.cm = DBH.23*2.54)%>%
  group_by(BLOCK, PLOT, CODE) %>%
  summarize(
    bapa = sum(ba * ef),                           
    tpa = sum(ef),                                
    qmd = qmd(bapa, tpa),                          
    rd = relative.density.index(bapa, qmd),        
    volume.ha = sum(new.vol.m3*ef*2.47, na.rm = TRUE),        
    CCF = (sum(c.area*ef, na.rm = TRUE)/43560)*100,                 
    LAI = mean(LAI, na.rm = TRUE))

head(plot.summary)

bl.t <- plot.summary%>%
  left_join(.,site)%>%
  mutate(baph = bapa/4.356,
         tph = tpa*2.47,
         qmd = qmd*2.54)
head(bl.t)

names(bl.t)

require(MASS)
b <- boxcox(lm(volume.ha ~ 1,data=bl.t))
# Exact lambda
lambda <- b$x[which.max(b$y)]
lambda #lambda is now 1.353535 so box cox transformation is no longer needed

library(fitdistrplus)
library(agricolae)
require(leaps)

descdist(bl.t$volume.ha, boot = 1000)
plot(density(bl.t$volume.ha))
plot(density(exp(bl.t$volume.ha)))
plot(density(sqrt(bl.t$volume.ha)))
plot(density(log(bl.t$volume.ha)))

#bl.t$sq.vol <- sqrt(bl.t$volume.ha)

names(bl.t)
modz <- regsubsets(volume.ha~baph+tph+qmd+rd+CCF+LAI+BS_Suitability+WS_Suitability+RS_Suitability+
                     elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+vpdmin+McNab+Bolstad+
                     Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+WD2000+WHC+ex.mg+ex.ca+ex.k+ph+dep+nit+SWC2+
                     NS.rd+WS.rd+BS.rd+RS.rd+NS.vol+WS.vol+RS.vol+BS.vol,
                   data=bl.t)
summary(modz)

modz2 <- regsubsets(volume.ha~LAI+BS_Suitability+WS_Suitability+RS_Suitability+
                      elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+vpdmin+McNab+Bolstad+
                      Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+WD2000+WHC+ex.mg+ex.ca+ex.k+ph+dep+nit+SWC2+
                      NS.rd+WS.rd+BS.rd+RS.rd+NS.vol+WS.vol+RS.vol+BS.vol,
                    data=bl.t)
summary(modz2)

#bapa + tpa + qmd + rd + CCF

ocean <- lm(volume.ha~qmd+rd+McNab+Winds10+NS.rd+BS.rd, data=bl.t)
summary(ocean)
vif(ocean)
AIC(ocean)

lake <- lm(volume.ha~rd+dep*CODE, data=bl.t)
summary(lake)  
vif(lake)
AIC(lake)

pond <- lme(volume.ha~rd+NS.rd+dep, 
            random=~1|BLOCK,
            data=bl.t)

summary(pond)
performance(pond)


testmod<-lme(volume.ha~CODE,
             random=~1|BLOCK,
             data=bl.t)

summary(testmod)
performance(testmod)

river <- lm(volume.ha~rd+NS.rd+dep, data = bl.t)
summary(river)
vif(river)
performance(river)

MB <- mean(residuals(pond))
print(paste("Mean Bias (MB):", MB))

MAB <- mean(abs(residuals(pond)))
print(paste("Mean Absolute Bias (MAB):", MAB))

bl.t$fit.vol <- predict(pond, newdata = bl.t)
plot(bl.t$fit.vol,bl.t$volume.ha)

f1 <- lme(fit.vol ~ CODE, random = ~1 | BLOCK, data = bl.t)
f1

library(emmeans)
emm <- emmeans(f1, ~ CODE)
pairs(emm, adjust = "tukey")

library(multcomp)
library(multcompView)
cld(emm, Letters = letters)


bl.t %>%
  group_by(CODE) %>%
  summarise(
    mean_fit_vol = mean(fit.vol, na.rm = TRUE),
    sd_fit_vol = sd(fit.vol, na.rm = TRUE),
    min_fit_vol = min(fit.vol, na.rm = TRUE),
    max_fit_vol = max(fit.vol, na.rm = TRUE)
  )


#png("~/Desktop/SMC_Sinuosity_Model_Output.png",units='in',height=5.5,width=14,res=1000)
#theme_set(theme_bw(16))

#library(ggplot2)
#ggplot(mydf2,aes(x=x,y=predicted,colour=group))+
  #geom_line(aes(linetype=group,color=group),size=1)+
  #labs(x="Tree height (m)",y="Probabilty of sinuosity (%)")+
  #ylim(0,55)+
  #xlim(4,6)+
  #facet_wrap(~facet)+
  #theme_bw(18) 
  #theme(legend.position="none")
  #scale_y_continuous(trans="sq")
  #scale_color_manual(values=c('gray0','gray70','gray40'))+
  #scale_fill_manual(values=c('gray0','gray70','gray40'), name="fill")

#back transform volume
#mydf2 <- ggpredict(full.m1, terms = c("rd", "CODE","aspect"))

# Back-transform from square root
#mydf2$predicted <- mydf2$predicted^2
#mydf2$conf.low <- mydf2$conf.low^2
#mydf2$conf.high <- mydf2$conf.high^2

#ggplot(mydf2, aes(x = x, y = predicted, colour = group)) +
  #geom_line(aes(linetype = group, color = group), size = 1) +
  #facet_wrap(~facet) +
  #ggtitle(expression(paste("Aspect (", degree, ")"))) +
  #labs(
    #x = "Relative Density (%)",
    #y = expression(paste("Volume (", m^3, " ", ha^-1, ")")),
    #colour = "Treatment",
    #linetype = "Treatment"
  #) +
  #coord_cartesian(ylim = c(0, 300)) +
  #theme_bw(18) +
  #theme(
    #plot.title = element_text(hjust = 0.5)  # centers the title
  #)


# Back-transform volume
#mydf2 <- ggpredict(full.m1, terms = c("rd", "CODE"))

# Back-transform from square root
#mydf2$predicted <- mydf2$predicted^2
#mydf2$conf.low <- mydf2$conf.low^2
#mydf2$conf.high <- mydf2$conf.high^2

# Rename columns for clarity (if necessary)
# str(mydf2) to check names; assuming 'x' = rd and 'group' = CODE
# If not, adjust accordingly

#ggplot(mydf2, aes(x = x, y = predicted, colour = group)) +
  #geom_line(aes(linetype = group), size = 1) +
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, colour = NA) +
  #ggtitle(expression(paste("Aspect (", degree, ")"))) +
  #labs(
    #x = "Relative Density (%)",
    #y = expression(paste("Volume (", m^3, " ", ha^{-1}, ")")),
    #colour = "Treatment",
    #linetype = "Treatment",
    #fill = "Treatment"
  #) +
  #coord_cartesian(ylim = c(0, 300)) +
  #theme_bw(base_size = 18) +
  #theme(
    #plot.title = element_text(hjust = 0.5)
  #)

#-------------------------------------------------------------------------------
# OY & TOY
#-------------------------------------------------------------------------------
# Overyielding
# Volumes for mixture & monoculture counterparts
mix <- c(120,93,73)
mono1 <- c(107,89,106)
mono2 <- c(79,84,82)

mean.mix <- mean(mix)
mean.mono1 <- mean(mono1)
mean.mono2 <- mean(mono2)
mean.mono.avg <- mean((mono1 + mono2) / 2)

sd.mix <- sd(mix)
sd.mono.avg <- sd((mono1 + mono2) / 2)

overyielding <- mean.mix / mean.mono.avg

print(overyielding)
print(sd.mix)
print(sd.mono.avg)

# combine data for ANOVA
data <- c(mix, (mono1 + mono2) / 2)
group <- factor(rep(c("mix", "mono avg"), each = length(mix)))

# perform 1-way ANOVA
anova_result <- aov(data ~ group)
summary(anova_result)

# Transgressive Overyielding

# Mixture (e.g., BR) and best monoculture (e.g., B)
mix <- c(70,67,71)            
bestmono <- c(82,84,79)

mean.mix <- mean(mix)
mean.bestmono <- mean(bestmono)

sd.mix <- sd(mix)
sd.bestmono <- sd(bestmono)

toy <- mean.mix / mean.bestmono

print(mean.mix)
print(mean.bestmono)
print(sd.mix)
print(sd.mono)
print(toy)

# combine data
data <- c(mix, bestmono)
group <- factor(c(rep("Mixture", length(mix)), rep("BestMono", length(bestmono))))

# perform 1-way ANOVA
anova_result <- aov(data ~ group)
summary(anova_result)

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
  filter(!(CODE == "BW" & SPP == "RS"),   
         !(CODE == "NR" & SPP == "BS"),   
         !(CODE == "RW" & SPP %in% c("NS", "BS")))  

ggplot(spruce.dbh2, aes(x = DBH.23, fill = SPP, color = SPP)) +
  geom_density(alpha = 0.5) +  
  labs(x = "DBH (in)", y = "Density") +
  scale_x_continuous(breaks = seq(0, max(spruce2$DBH.23, na.rm = TRUE), by = 2),
                     expand = expansion(mult = c(0.05, 0.1))) +  
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~ CODE)

ggplot(spruce.dbh2, aes(x = DBH.23, fill = SPP, color = SPP)) +
  geom_density(aes(y = ..count..), adjust = 1, alpha = 0.3) + 
  labs(x = "DBH (in)", y = "Trees per Acre", fill = "Species", color = "Species") + 
  scale_x_continuous(breaks = seq(0, max(spruce.filtered$DBH.23, na.rm = TRUE), by = 2),
                     expand = expansion(mult = c(0.05, 0.1))) +  
  theme_minimal() +
  theme(legend.position = "right",  
        legend.title = element_text(size = 10)) +  
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
spruce.dbh3 <- spruce.dbh2 %>% 
  filter(CODE %in% c("BN", "BR", "BW", "NR", "NW", "RW"))

ggplot(spruce.dbh3, aes(x = DBH.23 * 2.54, fill = SPP, color = SPP)) +
  geom_density(alpha = 0.5) +  
  labs(x = "DBH (cm)", y = "Density") +  
  theme_minimal() +
 theme(legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~ CODE)

#DBHxTPA
spruce.dbh2 %>%
  mutate(DBH_cm = DBH * 2.54,  
         DBH_class = cut(DBH_cm, breaks = seq(0, max(DBH_cm), by = 3), include.lowest = TRUE)) %>%  
  filter(!is.na(DBH_class)) %>%  
  group_by(CODE, DBH_class) %>%
  summarise(count = n()) %>%
  mutate(count = count * 10 * 2.47105) %>%  
  ggplot(aes(x = DBH_class, y = count)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ CODE) +
  labs(x = "DBH (cm)", y = "Trees per Hectare") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

spruce.dbh2 %>%
  filter(!CODE %in% c("B", "N", "W", "R")) %>%  
  filter(!(CODE == "BN" & SPP == "WS")) %>%  
  mutate(DBH_cm = DBH * 2.54,  
         DBH_class = cut(DBH_cm, breaks = seq(0, max(DBH_cm), by = 3), include.lowest = TRUE)) %>%
  filter(!is.na(DBH_class)) %>% 
  group_by(CODE, SPP, DBH_class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(count = count * 10 * 2.47105) %>% #10 for exp factor, 2.47 from ac to ha
  ggplot(aes(x = DBH_class, y = count, fill = SPP)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ CODE) +
  labs(x = "DBH (cm)", y = "Trees per Hectare") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")

spruce.dbh2 %>%
  filter(!CODE %in% c("B", "N", "W", "R")) %>%  
  filter(!(CODE == "BN" & SPP == "WS")) %>%  
  mutate(DBH_cm = DBH * 2.54) %>%
  filter(!is.na(DBH_cm)) %>%
  ggplot(aes(x = DBH_cm, fill = SPP)) +
  geom_histogram(binwidth = 3, aes(y = ..count.. * 10 * 2.47105), alpha = 0.4, position = "identity") +  
  facet_wrap(~ CODE) +
  labs(x = "DBH (cm)", y = "Trees per Hectare") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")

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

ss.df <- ss.df %>% filter(CODE != "C")

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
#NR aov shape p-value: 0.00334 *
#   aov scale p-value: 0.000535 *
#   tukey shape: a b
#   tukey scale: a b
#NW aov shape p-value: 0.0815
#   aov scale p-value: 0.367
#   tukey shape: a a
#   tukey scale: a a
#RW aov shape p-value: 0.0155
#   aov scale p-value: 4e-04 *
#   tukey shape: a a
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

library(multcomp)

ss.df$SPP <- as.factor(ss.df$SPP)
ss.df$CODE <- as.factor(ss.df$CODE)

amod <- lm(Shape ~ SPP + CODE, data = ss.df)  
bmod <- lm(Scale ~ SPP + CODE, data = ss.df)  

summary(amod)
summary(bmod)

shape.glht <- glht(amod, linfct = mcp(CODE = "Tukey")) #just accounts for differences amongst CODEs 
summary(shape.glht)

scale.glht <- glht(bmod, linfct = mcp(CODE = "Tukey"))
summary(scale.glht)

# test for differences amongst SPP in each CODE individually

#Shape
amod <- lm(Shape ~ SPP, data = ss.df[ss.df$CODE == "RW", ])  
shape.glht <- glht(amod, linfct = mcp(SPP = "Tukey"))
summary(shape.glht)

#Scale
bmod <- lm(Scale ~ SPP, data = ss.df[ss.df$CODE == "RW", ])  
scale.glht <- glht(bmod, linfct = mcp(SPP = "Tukey"))
summary(scale.glht)


# test for differences amongst SPP in each CODE all at once

ss <- ss.df %>% filter(!CODE %in% c("B", "N", "R", "W"))

amod <- lm(Shape ~ SPP * CODE, data = ss)  
bmod <- lm(Scale ~ SPP * CODE, data = ss)  

unique.code <- unique(ss$CODE)

#loop isnt working correctly
for (code in unique.code) {
  amod <- lm(Shape ~ SPP, data = ss[ss$CODE == code, ]) 
  shape.glht <- glht(amod, linfct = mcp(SPP = "Tukey"))
  cat("\n\nTukey's HSD for Shape in CODE:", code, "\n")
  print(summary(shape_glht))
  
  bmod <- lm(Scale ~ SPP, data = ss[ss$CODE == code, ])
  scale.glht <- glht(bmod, linfct = mcp(SPP = "Tukey"))
  cat("\n\nTukey's HSD for Scale in CODE:", code, "\n")
  print(summary(scale_glht))
}


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


library(multcomp)

ss.ht.df$SPP <- as.factor(ss.ht.df$SPP)
ss.ht.df$CODE <- as.factor(ss.ht.df$CODE)

cmod <- lm(Shape ~ SPP + CODE, data = ss.ht.df)  
dmod <- lm(Scale ~ SPP + CODE, data = ss.ht.df)  

summary(cmod)
summary(dmod)

shape.glht <- glht(cmod, linfct = mcp(CODE = "Tukey")) #just accounts for differences amongst CODEs 
summary(shape.glht)

scale.glht <- glht(dmod, linfct = mcp(CODE = "Tukey"))
summary(scale.glht)

# test for differences amongst SPP in each CODE individually

#Shape
cmod <- lm(Shape ~ SPP, data = ss.ht.df[ss.ht.df$CODE == "RW", ])  
shape.glht <- glht(cmod, linfct = mcp(SPP = "Tukey"))
summary(shape.glht)

#Scale
dmod <- lm(Scale ~ SPP, data = ss.ht.df[ss.ht.df$CODE == "RW", ])  
scale.glht <- glht(dmod, linfct = mcp(SPP = "Tukey"))
summary(scale.glht)

#-------------------------------------------------------------------------------
# weevil damage
#-------------------------------------------------------------------------------
# filter for just NS stems and CODEs N, BN, NR, NW
christmas <- picea %>%
  filter(SPP == "NS", CODE %in% c("N", "NW", "NR", "BN"))

# eeshape DAM1, DAM2, DAM3 into long format
newyears <- christmas %>%
  pivot_longer(cols = c(DAM1, DAM2, DAM3), names_to = "DAM_num", values_to = "DAM_value") %>%
  filter(DAM_value %in% c("FK/WV", "WV"))

# Each row now represents one tree with DAM_value == FK/WV or WV
# Count by BLOCK, PLOT, and CODE
weevil<- newyears %>%
  group_by(CODE, BLOCK, PLOT, SPP) %>%
  summarise(n_trees = n(), .groups = "drop")

w1 <- lm(n_trees ~ CODE, data = weevil)
summary(w1)

library(emmeans)
emm <- emmeans(w1, ~ CODE)
pairs(emm, adjust = "tukey")

library(multcomp)
library(multcompView)
cld(emm, Letters = letters)



