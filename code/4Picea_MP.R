#BobSpruce Data Analysis
#Ashley Carter
#Forest Management Lab
# Edits by Premer, M.I. 20230907

#loadindata
BobSpruce <- read.csv('~/Desktop/BobSpruce.csv')
head(BobSpruce)


#loadpackages
library(dplyr)
library(mosaic)
library(forcats)
library(tidyverse)

BobSpruce$DBH <- as.numeric(BobSpruce$DBH)
BobSpruce$Ht <- as.numeric(BobSpruce$Ht)
BobSpruce[is.na(BobSpruce)] <- 0
plot(density(BobSpruce$Ht))
# something isn't right... a 300 ft tree? 
summary(BobSpruce)
BobSpruce$Ht <- ifelse(BobSpruce$Ht<200,BobSpruce$Ht,0)
library(lattice)
xyplot(Ht~DBH|Spp.,data=BobSpruce,na.omit=TRUE)

# there is no way that is ft... must have been collected in meters for height. 
Bob <- BobSpruce %>% 
  mutate(Ht.ft = Ht*3.28)%>%
  mutate(DBH.in = DBH*0.3937)%>%
  mutate(tf = 10)%>%
  mutate(ba = DBH.in^2*0.005454)%>%
  mutate(batf = tf*ba)

Rob <- Bob%>%
  group_by(Block,Plot)%>%
  summarize(bapa = sum(batf),
            tpa = sum(tf))

Robert <- Bob%>%
  group_by(Block,Plot,Spp.)%>%
  summarize(sp.bapa = sum(batf),
            sp.tpa = sum(tf))%>%
  ungroup()%>%
  left_join(.,Rob)%>%
  mutate(prop.tpa = sp.tpa/tpa)%>%
  mutate(prop.ba = sp.bapa/bapa)%>%
  mutate(IV = (prop.tpa+prop.ba)/2)
  
# Double check, all plots should sum to 1
check <- Robert%>%
  group_by(Block,Plot)%>%
  summarize(checker = sum(IV))
print(check,n=30)
# bueno..
print(Rob)

xyplot(bapa~Plot|Block*Spp.,data=Robert)
# end 20230907
##########################################
##########################################


#remove NA values in DBH
BobSpruce<-BobSpruce%>%
  drop_na(DBH)

#convert DBH column from character to numeric
BobSpruce.DBH <- as.numeric(as.character(BobSpruce$DBH))

#remove NA values in DBH
BobSpruce<-BobSpruce%>%
  drop_na(DBH)

#remove NA values in Ht
BobSpruce<-BobSpruce%>%
  drop_na(Ht)

#convert Ht column from character to numeric
BobSpruce.Ht <- as.numeric(as.character(BobSpruce$Ht))

#add in expansion factor
BobSpruce <- BobSpruce%>%
  mutate(PEF = 10)

#calculating basal area
BobSpruce <- BobSpruce%>%
  mutate(BA = (DBH^2)*0.005454)

#remove NA values in basal area
BobSpruce<-BobSpruce%>%
  drop_na(BA)

#basal area per acre by tree
BobSpruce<-BobSpruce%>%
  mutate(TREE_BA_AC = BA * PEF)

#basal area and tpa calculations for plot level by block
plot.summary <- BobSpruce%>%
  group_by(Block, Plot)%>%
  summarise(TPA_total = sum(PEF),
            bapa_total = sum(TREE_BA_AC))
view(plot.summary)

#basal area and tpa calculations by species
species.summary <- BobSpruce%>%
  group_by(Block, Plot, Spp.)%>%
  summarise(TPA = sum(PEF),
            bapa = sum(TREE_BA_AC))
view(species.summary)

#join species and plot tables together
BobSpruce <- left_join(species.summary, plot.summary)

#importance values
BobSpruce.prop <- BobSpruce%>%
  mutate(prop_tpa = (TPA/TPA_total),
         prop_ba = (bapa/bapa_total))
view(BobSpruce.prop)

BobSpruce.IV<- BobSpruce%>%
  mutate(iv = (prop_tpa+prop_ba)/2)
view(BobSpruce.IV)
