#install.packages("devtools")
#install.packages("G:/My Drive/MEForLab",
                 #repos = NULL,
                 #type = "source")
devtools::install_github("piketremor/MEForLab")
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

############### calculate species IVs
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

#view(picea.species)

############### Double check, all plots should sum to 1
check <- picea.species%>%
  group_by(uid)%>%
  summarize(checker = sum(IV))
print(check,n=30)
# they sum to 1, good

############### imputing tree heights
library(lme4)

#select for SITEid, PLOTid, and TREE to be in tree locations dataset
tree_species <- picea%>%
  select(Block, Plot, Tree, Species)

#filter for the living trees
picea <- filter(picea, StatusCode.2023 == 1)

#filter for heights above 0
height.frame <- filter(picea, HT.2023 > 0)

#create model, Henrikson Equation (be wary of tree > 5" DBH)
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
  mapply(vol.calc,SPP=tree_predict$Species,DBH=tree_predict$DBH.2023,HT=tree_predict$fin.ht)

xyplot(vol~DBH|SPP,data=tree_predict)
xyplot(fin.ht~DBH|SPP,data=tree_predict)
xplot(fin.ht~DBH,data=tree_predict)

############### basal area larger (bal)
spruce.bal <- picea%>%
  mutate(basal.area = picea$DBH.2023^2*0.005454)%>%
  mutate(bapa = basal.area*10)%>%
  group_by(uid)%>%
  arrange(desc(DBH.2023),.by_group = TRUE)%>%
  mutate(bal = lag(cumsum(bapa)))
spruce.bal$bal[is.na(spruce.bal$bal)] <- 0
xyplot(bal~DBH.2023|Species,data=spruce.bal)

################ height larger (htl)
spruce.htl <- tree_predict%>%
group_by(uid)%>%
  arrange(desc(fin.ht),.by_group = TRUE)%>%
  mutate(htl = lag(cumsum(fin.ht)))
spruce.htl$htl[is.na(spruce.htl$htl)] <- 0

xyplot(htl~fin.ht|Species,data=spruce)

############### ht/diameter ratios
ht_dbh <- tree_predict %>% 
  mutate(ht_dbh = fin.ht/DBH.2023)

############### vicary SI, error need to enter OS and NS
tree_predict["SI"] <- 
  mapply(vicary.site,SPP=tree_predict$Species,ht=tree_predict$fin.ht,age=28)


################ ANOVA: does diameter vary by species
aov.1<- aov(DBH.2023~Species, data=picea)
summary(aov.1)

TukeyHSD(aov.1)

#Homogeneity of variances
plot(aov.1, 1)

#Normality
plot(aov.1, 2)

############### ANOVA: does diameter vary by species grouped by id
picea <- picea%>%
  unite(ID, 
        Block, 
        Plot,
        remove = FALSE,
        sep = ".")

aov.2 <- aov(DBH.2023~Species+ID,data=picea)
summary(aov.2)

# if you were look for an interaction use formula below
# aov.3 <- aov(DBH.2023) ~ Species * id, data = my_data)

group_by(picea, Species, ID) %>%
  summarise(
    count = n(),
    mean = mean(DBH.2023, na.rm = TRUE),
    sd = sd(DBH.2023, na.rm = TRUE)
  )

# pairwise comparisons between groups
TukeyHSD(aov.2, which = "Species")
TukeyHSD(aov.2, which = "ID")

#Homogeneity of variances
plot(aov.2, 1)

library(car) #if p-value is <.05 then var between groups is significantly different, don't want that
leveneTest(DBH.2023 ~ Species*id, data = picea)

#Normality
plot(aov.2, 2)

############### ANOVA: does diameter vary by species mixture

############### ANOVA: does diameter vary by species mixture grouped by id

