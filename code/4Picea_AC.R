#install.packages("devtools")
#install.packages("G:/My Drive/MEForLab",
                 #repos = NULL,
                 #type = "source")

library(MEForLab)
library(devtools)
library(dplyr)
library(mosaic)
library(forcats)
library(tidyverse)

setwd("G:/Shared drives/4Picea/4Picea/raw")
picea <- read.csv("4Picea.Data.csv")

picea$DBH.2023 <- as.numeric(picea$DBH.2023)
picea$DBH.2023[is.na(picea$DBH.2023)] <- 0
picea$uid <- paste0(picea$Block,".",picea$Plot)

picea$tree.factor <- 10

picea.summary <- picea%>%
  filter(.,StatusCode.2023=="1")%>%
  mutate(tree.ba = DBH.2023^2*0.005454)%>%
  group_by(uid)%>%
  summarize(bapa = sum(tree.ba*tree.factor), # tree.factor is your expansion factor
            tpa = sum(tree.factor))

plot(picea.summary$bapa,picea.summary$tpa)
# most plots are between 20 and 140 bapa

# calculate species IVs
picea.species <- picea%>%
  filter(.,StatusCode.2023=="1")%>%
  group_by(uid,Species)%>%
  summarize(sp.bapa = sum(DBH.2023^2*.005454*tree.factor),
            sp.tpa = sum(tree.factor))%>%
  ungroup()%>%
  left_join(.,picea.summary)%>%
  mutate(prop.tpa = sp.tpa/tpa)%>%
  mutate(prop.ba = sp.bapa/bapa)%>%
  mutate(IV = (prop.tpa+prop.ba)/2)

# Double check, all plots should sum to 1
check <- picea.species%>%
  group_by(uid)%>%
  summarize(checker = sum(IV))
print(check,n=30)
# they sum to 1, good

#############################################

#imputing tree heights
library(lme4)

#select for SITEid, PLOTid, and TREE to be in tree locations dataset
tree_species <- picea%>%
  select(Block, Plot, Tree, Species)

#filter for the living trees
picea <- filter(picea, StatusCode.2023 == 1)

#filter for heights above 0
height.frame <- filter(picea, HT.2023 > 0)

#create model
ht.lm <- lmer(HT.2023 ~ log(DBH.2023) + (1 | Species/Plot/Block), data = height.frame)
summary(ht.lm)

#predict heights
tree_predict <- picea %>% 
  mutate(PRD_HT = predict(ht.lm, picea, re.form = NULL, allow.new.levels = TRUE))

tree_predict$HT.2023[is.na(tree_predict$HT.2023)] <- 0
tree_predict$PRD_HT[is.na(tree_predict$PRD_HT)] <- 0
tree_predict$fin.ht <- ifelse(tree_predict$HT.2023<1,tree_predict$PRD_HT,tree_predict$HT.2023)
tree_predict$fin.ht <- ifelse(tree_predict$fin.ht<0,0,tree_predict$fin.ht)

tree_predict$Species[is.na(tree_predict$Species)] <- "OS"

### getting an error here, not sure why
tree_predict["vol"] <- 
  mapply(vol_calc,SPP=tree_predict$Species,DBH=tree_predict$DBH.2023,HT=tree_predict$fin.ht)


xyplot(vol~DBH|SPP,data=tree_predict)

xyplot(fin.ht~DBH|SPP,data=tree_predict)

xplot(fin.ht~DBH,data=tree_predict)

########################################

#ht/diameter ratios

ht_dbh <- tree_predict %>% 
  mutate(ht_dbh = fin.ht/DBH.2023)
 
#next summarize by block, plot, species?
   
