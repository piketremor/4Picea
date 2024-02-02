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
picea <- read.csv("4Picea.csv")

picea$DBH.23 <- as.numeric(picea$DBH.23)
picea$DBH.23[is.na(picea$DBH.23)] <- 0
picea$uid <- paste0(picea$BLOCK,".",picea$PLOT)

picea$tree.factor <- 10

picea.summary <- picea%>%
  filter(.,STATUS.23=="1")%>%
  mutate(tree.ba = DBH.23^2*0.005454)%>%
  group_by(uid)%>%
  summarize(bapa = sum(tree.ba*tree.factor), # tree.factor is your expansion factor
            tpa = sum(tree.factor))

plot(picea.summary$bapa,picea.summary$tpa) # most plots are between 20 and 140 bapa

######################### calculate species IVs

picea.species <- picea%>%
  filter(.,STATUS.23=="1")%>%
  group_by(uid,SPP)%>%
  summarize(sp.bapa = sum(DBH.23^2*.005454*tree.factor),
            sp.tpa = sum(tree.factor))%>%
  ungroup()%>%
  left_join(.,picea.summary)%>%
  mutate(prop.tpa = sp.tpa/tpa)%>%
  mutate(prop.ba = sp.bapa/bapa)%>%
  mutate(IV = (prop.tpa+prop.ba)/2)

#view(picea.species)

####################### double check, all plots should sum to 1
check <- picea.species%>%
  group_by(uid)%>%
  summarize(checker = sum(IV))
print(check,n=30) # they sum to 1, good

###################### imputing tree heights
library(lme4)

# select for SITEid, PLOTid, and TREE to be in tree locations dataset
tree_species <- picea%>%
  select(BLOCK, PLOT, TREE, SPP)

# filter for the living trees
picea <- filter(picea, STATUS.23 == 1)

# filter for heights above 0
height.frame <- filter(picea, HT.23 > 0)

# create model, Henrikson Equation (be wary of tree > 5" DBH)
ht.lm <- lmer(HT.23 ~ log(DBH.23) + (1 | SPP/PLOT/BLOCK), data = height.frame) #error means that model is likely overfitted
summary(ht.lm)

# predict heights
tree_predict <- picea %>% 
  mutate(PRD_HT = predict(ht.lm, picea, re.form = NULL, allow.new.levels = TRUE))

tree_predict$HT.23[is.na(tree_predict$HT.23)] <- 0
tree_predict$PRD_HT[is.na(tree_predict$PRD_HT)] <- 0
tree_predict$fin.ht <- ifelse(tree_predict$HT.23<1,tree_predict$PRD_HT,tree_predict$HT.23)
tree_predict$fin.ht <- ifelse(tree_predict$fin.ht<0,0,tree_predict$fin.ht)

tree_predict$SPP[is.na(tree_predict$SPP)] <- "OS"

## getting an error here, not sure why
tree_predict["vol"] <- 
  mapply(vol.calc,SPP=tree_predict$SPP,DBH=tree_predict$DBH.23,HT=tree_predict$fin.ht)

xyplot(vol~DBH|SPP,data=tree_predict)
xyplot(fin.ht~DBH|SPP,data=tree_predict)
xplot(fin.ht~DBH,data=tree_predict)

######################## basal area larger (bal)
spruce.bal <- picea%>%
  mutate(basal.area = picea$DBH.23^2*0.005454)%>%
  mutate(bapa = basal.area*10)%>%
  group_by(uid)%>%
  arrange(desc(DBH.23),.by_group = TRUE)%>%
  mutate(bal = lag(cumsum(bapa)))
spruce.bal$bal[is.na(spruce.bal$bal)] <- 0
xyplot(bal~DBH.23|SPP,data=spruce.bal)

###################### height larger (htl)
spruce.htl <- tree_predict%>%
group_by(uid)%>%
  arrange(desc(fin.ht),.by_group = TRUE)%>%
  mutate(htl = lag(cumsum(fin.ht)))
spruce.htl$htl[is.na(spruce.htl$htl)] <- 0

xyplot(htl~fin.ht|SPP,data=spruce.htl)

###################### ht/diameter ratios
ht_dbh <- tree_predict %>% 
  mutate(ht_dbh = fin.ht/DBH.23)

###################### vicary SI, error need to enter OS and NS
#tree_predict["SI"] <- 
  #mapply(vicary.site,SPP=tree_predict$Species,ht=tree_predict$fin.ht,age=28)


##################### ANOVA: does diameter vary by species (one-way)
aov.1<- aov(DBH.23~SPP, data=picea)
summary(aov.1)

TukeyHSD(aov.1)

#Homogeneity of variances
plot(aov.1, 1)

#Normality
plot(aov.1, 2)

############### ANOVA: does diameter vary by species grouped by id (two-way)
picea <- picea%>%
  unite(ID, 
        BLOCK, 
        PLOT,
        remove = FALSE,
        sep = ".")

aov.2 <- aov(DBH.23~SPP+ID,data=picea)
summary(aov.2)

# if you were looking for an interaction use formula below
#aov.3 <- aov(DBH.2023) ~ Species * id, data = my_data)

group_by(picea, SPP, ID) %>%
  summarise(
    count = n(),
    mean = mean(DBH.23, na.rm = TRUE),
    sd = sd(DBH.23, na.rm = TRUE)
  )

# pairwise comparisons between groups
TukeyHSD(aov.2, which = "SPP")
TukeyHSD(aov.2, which = "ID")

#Homogeneity of variances
plot(aov.2, 1)

library(car) #if p-value is <.05 then var between groups is significantly different, don't want that
leveneTest(DBH.23 ~ SPP*ID, data = picea)

#Normality
plot(aov.2, 2)

##################### ANOVA: does diameter vary by pure vs mixed (one-way)
aov.3<- aov(DBH.23~CODE, data=picea)
summary(aov.3)

TukeyHSD(aov.3)

#Homogeneity of variances
plot(aov.3, 1)

#Normality
plot(aov.3, 2)

################# ANOVA: does diameter vary by mixed vs pure grouped by ID (two-way)
picea <- picea%>%
  unite(ID, 
        BLOCK, 
        PLOT,
        remove = FALSE,
        sep = ".")

aov.4 <- aov(DBH.23~CODE+ID,data=picea)
summary(aov.4)

# if you were looking for an interaction use formula below
#aov.3 <- aov(DBH.2023) ~ CODE * id, data = my_data)

group_by(picea, CODE, ID) %>%
  summarise(
    count = n(),
    mean = mean(DBH.23, na.rm = TRUE),
    sd = sd(DBH.23, na.rm = TRUE)
  )

# pairwise comparisons between groups
TukeyHSD(aov.4, which = "CODE")
TukeyHSD(aov.4, which = "ID")

#Homogeneity of variances
plot(aov.4, 1)

library(car) #if p-value is <.05 then var between groups is significantly different, don't want that
leveneTest(DBH.23 ~ CODE*ID, data = picea)

#Normality
plot(aov.4, 2)

############## Matt Russell ANOVA Textbook Exercises
p.diameter <- ggplot(picea, aes(factor(SPP), DBH.23)) + 
  geom_boxplot()+
  ylab("Diameter (inches)") +
  xlab("Species")
p.diameter

# one way diameter + species
#picea <- picea %>% 
#mutate(Species.fact = as.factor(Species))

picea.aov <- lm(DBH.23 ~ SPP, data = picea)

anova(picea.aov)

#library(broom)
#library(patchwork)

#picea.resid <- diagPlot(picea.aov) #no diagPlot function anymore?

#p.diag <-delta

#(p.diag$p.resid | p.diag$p.stdresid) /
#(p.diag$p.srstdresid | p.diag$p.cooks)

pairwise.t.test(picea$DBH.23, picea$SPP, p.adj = "bonferroni")

library(agricolae)
lsd.picea <- LSD.test(picea.aov, "SPP", p.adj = "bonferroni")
lsd.picea$groups

picea_summ <- picea %>% 
  group_by(SPP)  %>%  
  summarize(n.diameter = n(),
            mean.diameter = mean(DBH.23),
            sd.diameter = sd(DBH.23))
picea_summ

limits <- aes(ymax = mean.diameter + sd.diameter, 
              ymin = mean.diameter - sd.diameter)
p.diameter <- ggplot(picea_summ, aes(SPP, mean.diameter)) +
  geom_bar(stat = "identity") +
  geom_errorbar(limits, width = 0.25) +
  ylab("Diameter (Inches") +
  xlab("Species") +
  geom_text(aes(label = c("a", "b", "c", "d", "d", "d")), vjust = -6) +
  scale_y_continuous(limits = c(0, 0.3))
p.diameter #plot doesn't work, not sure why

# one way ANOVA diameter x species mix
p.diameter <- ggplot(picea, aes(factor(CODE), DBH.23)) + 
  geom_boxplot()+
  ylab("Diameter (inches)") +
  xlab("Species Mix")
p.diameter
piceamix.aov <- lm(DBH.23 ~ CODE, data = picea)

anova(piceamix.aov)

pairwise.t.test(picea$DBH.23, picea$CODE, p.adj = "bonferroni")

library(agricolae)
lsd.picea <- LSD.test(piceamix.aov, "CODE", p.adj = "bonferroni")
lsd.picea$groups

picea_summ2 <- picea %>% 
  group_by(CODE)  %>%  
  summarize(n.diameter = n(),
            mean.diameter = mean(DBH.23),
            sd.diameter = sd(DBH.23))
picea_summ2

# two way diameter x species + ID
ggplot(picea, aes(SPP, DBH.23, fill = ID)) +
  geom_boxplot() +
  ylab("Diameter (inches)") +
  xlab("Species")

picea.aov2 <- lm(DBH.23 ~ SPP * ID, data = picea)
anova(picea.aov2)

pairwise.t.test(picea$DBH.23, picea$ID, p.adj = "bonferroni")

lsd.picea2 <- LSD.test(picea.aov2, "ID", p.adj = "bonferroni")
lsd.picea2$groups

picea_summ2 <- picea %>% 
  group_by(ID) %>%  
  summarize(n.diameter = n(),
            mean.diameter = mean(DBH.23),
            sd.diameter = sd(DBH.2023)) %>% 
  mutate(se.diameter = sd.diameter/sqrt(n.diameter))

picea_summ2

############## Matt Russell Linear Mixed Models Textbook Exercises

# generalized linear regression model for DBH + HT
ggplot(data=tree_predict, aes(x=DBH.23, y=fin.ht)) +
  geom_point() +
  geom_smooth(method = "lm")

cor <- cor.test(tree_predict$DBH.23, tree_predict$fin.ht, 
                method = "pearson")
cor #high correlation

lmodel <- lm(sqrt(DBH.23) ~ sqrt(fin.ht), data = tree_predict)
summary(lmodel)

AIC(lmodel)

# linear mixed models
library(nlme)
library(lme4)

# random effects on intercept
n_distinct(tree_predict$SPP)

ggplot(data=tree_predict, aes(x=DBH.23, y=fin.ht)) +
  geom_point() +
  facet_wrap(~SPP, ncol = 3) +
  labs(x = "Diameter at breast height (inches)",
       y = "Height (feet)")

picea.lme <- lmer(fin.ht ~ DBH.23 + (1 | SPP),
                    data = tree_predict)
summary(picea.lme)
ranef(picea.lme)

ggplot(tree_predict, aes(DBH.23, fin.ht)) +
  geom_point(size = 0.2) +
  geom_line(aes(y = predict(picea.lme), 
                group = SPP, 
                color = Species)) +
  labs(x = "Diameter at breast height (inches)",
       y = "Height (feet)")

coef(picea.lme)

#HT predictions with and without random effects
tree_predict <- tree_predict %>% 
  mutate(HT_pred_fixed_re = predict(picea.lme, 
                                    tree_predict, 
                                    re.form = NULL),
         HT_pred_fixed = predict(picea.lme, 
                                 tree_predict, 
                                 re.form = NA))

tree_predict %>% 
  top_n(10)

#random effects on slope, can introduce model complexity
picea.lme2 <- lmer(fin.ht ~ 1 + DBH.23 + (1 + DBH.23 | SPP),
                     data = tree_predict) 
summary(picea.lme2)

ggplot(tree_predict, aes(DBH.23, fin.ht)) +
  geom_point(size = 0.2) +
  geom_line(aes(y = predict(picea.lme2), 
                group = SPP, 
                color = SPP)) +
  labs(x = "Diameter at breast height (inches)",
       y = "Height (feet)")

# nested random effects on intercept
tree_predict <- tree_predict%>%
  unite(ID, 
        BLOCK, 
        PLOT,
        remove = FALSE,
        sep = ".")

n_distinct(tree_predict$ID)

tree_predict %>% #error found, doesn't work
  filter(ID == c(B3.2, B3.5)) %>% 
  ggplot(aes(DBH.23, fin.ht)) +
  geom_point() +
  facet_wrap(~ID) +
  labs(title = "Spruce at PEF",
       subtitle = "HT-DBH by plot ID",
       x = "Diameter at breast height (inches)",
       y = "Height (feet)")

picea.lme3 <- lmer(fin.ht ~ DBH.23 + (1 | SPP/ID),
                     data = tree_predict)
summary(picea.lme3)
ranef.lme3 <- ranef(picea.lme3)
plot(ranef.lme3)

AIC(lmodel, picea.lme, picea.lme3)









