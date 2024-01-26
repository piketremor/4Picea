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
#view(picea.species)

# Double check, all plots should sum to 1
check <- picea.species%>%
  group_by(uid)%>%
  summarize(checker = sum(IV))
print(check,n=30)
# they do not sum to 1, range from 1.01-1.18

#############################################
spruce <- picea%>%
  mutate(basal.area = picea$DBH.2023^2*0.005454)%>%
  mutate(bapa = basal.area*10)%>%
  group_by(uid)%>%
  arrange(desc(DBH.2023),.by_group = TRUE)%>%
  mutate(bal = lag(cumsum(bapa)))
spruce$bal[is.na(spruce$bal)] <- 0


