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

#-------------------------------------------------------------------------------
#clean data
#-------------------------------------------------------------------------------
picea$DBH.23 <- as.numeric(picea$DBH.23)
picea$DBH.23[is.na(picea$DBH.23)] <- 0
picea$uid <- paste0(picea$BLOCK,".",picea$PLOT)

#-------------------------------------------------------------------------------
#summary statistics
#-------------------------------------------------------------------------------
picea$tree.factor <- 10

picea <- picea%>%
  filter(.,STATUS.23=="1")%>%
  mutate(tree.ba = DBH.23^2*0.005454)%>%
  group_by(uid)%>%
  mutate(bapa = sum(tree.ba*tree.factor), # tree.factor is your expansion factor
            tpa = sum(tree.factor))

plot(picea.summary$bapa,picea.summary$tpa) # most plots are between 20 and 140 bapa
#-------------------------------------------------------------------------------
#calculate species IVs
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
#double check, all plots should sum to 1
#-------------------------------------------------------------------------------
check <- picea.species%>%
  group_by(uid)%>%
  summarize(checker = sum(IV))
print(check,n=30) # they sum to 1, good

#-------------------------------------------------------------------------------
#imputing tree heights (Henrikson equation)
#-------------------------------------------------------------------------------
#library(lme4)

# select for BLOCK, PLOT, TREE, SPP from picea dataset
#tree_species <- picea%>%
  #select(BLOCK, PLOT, TREE, SPP)

# filter for the living trees
#picea <- filter(picea, STATUS.23 == 1)

# filter for heights above 0
#height.frame <- filter(picea, HT.23 > 0)

# create model, Henrikson Equation (be wary of tree > 5" DBH)
#ht.lm <- lmer(HT.23 ~ log(DBH.23) + (1 | SPP/PLOT/BLOCK), data = height.frame) #error means that model is likely overfitted
#summary(ht.lm)

# predict heights
#tree_predict <- picea %>% 
  #mutate(PRD_HT = predict(ht.lm, picea, re.form = NULL, allow.new.levels = TRUE))

#tree_predict$HT.23[is.na(tree_predict$HT.23)] <- 0
t#ree_predict$PRD_HT[is.na(tree_predict$PRD_HT)] <- 0
#tree_predict$fin.ht <- ifelse(tree_predict$HT.23<1,tree_predict$PRD_HT,tree_predict$HT.23)
#tree_predict$fin.ht <- ifelse(tree_predict$fin.ht<0,0,tree_predict$fin.ht)

#tree_predict$SPP[is.na(tree_predict$SPP)] <- "OS"

#tree_predict["vol"] <- 
  #mapply(vol.calc,SPP=tree_predict$SPP,DBH=tree_predict$DBH.23,HT=tree_predict$fin.ht)

#xyplot(vol~DBH.23|SPP,data=tree_predict)
#xyplot(fin.ht~DBH.23|SPP,data=tree_predict)
#xplot(fin.ht~DBH.23,data=tree_predict)

#-------------------------------------------------------------------------------
#generating HT model
#-------------------------------------------------------------------------------
library(lme4)
mixer <- lmer(HT.23~log(DBH.23+0.1)+SPP+(1|PLOT/BLOCK), data=picea)
mixer2 <- lmer(HT.23~log(DBH.23+0.1)+(1|SPP/PLOT/BLOCK), data=picea)
AIC(mixer2, mixer)

library(nlme)
mixer3 <- lme(HT.23~log(DBH.23+0.1)+SPP, random=~1|PLOT/BLOCK, data=picea, na.action=na.omit)
AIC(mixer3)
summary(mixer3)
picea$BLOCK <- as.factor(picea$BLOCK)
picea$PLOT <- as.factor(picea$PLOT)

ht.mod2 <- nlme(HT.23~4.5+exp(a+(b/DBH.23+1)),
                data=picea,
                fixed=a+b~1,
                random=a+b~1|BLOCK/PLOT/SPP,
                na.action=na.pass,verbose=T,
                start=c(a=4.5, b=-6),
                control=nlmeControl(returnObject = TRUE,msMaxIter = 10000,maxIter = 5000))
summary(ht.mod2)

#-------------------------------------------------------------------------------
#imputing tree heights (Wykoff)
#-------------------------------------------------------------------------------
library(MEForLab)

picea$wykoff.ht <- predict(ht.mod2,picea)
picea$HT.23[is.na(picea$HT.23)] <- 0
picea$final.ht <- ifelse(picea$HT.23>0,picea$HT.23,picea$wykoff.ht)
xyplot(final.ht~DBH.23|SPP,data=picea)

#-------------------------------------------------------------------------------
#volume calculation
#-------------------------------------------------------------------------------
picea <- picea%>%
  mutate(vol=mapply(vol.calc,SPP=SPP,DBH=DBH.23,HT=final.ht))
xyplot(vol~DBH.23|SPP,data=picea)

#-------------------------------------------------------------------------------
#basal area larger (bal)
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

#-------------------------------------------------------------------------------
#height larger (htl)
#-------------------------------------------------------------------------------
picea <- picea%>%
group_by(uid)%>%
  arrange(desc(final.ht),.by_group = TRUE)%>%
  mutate(htl = lag(cumsum(final.ht)))
picea$htl[is.na(spruce.htl$htl)] <- 0
xyplot(htl~final.ht|SPP,data=spruce.htl)
xyplot(htl~final.ht|CODE,data=spruce.htl)

#-------------------------------------------------------------------------------
#height/diameter ratios
#-------------------------------------------------------------------------------
picea <- picea %>% 
  mutate(ht.dbh = final.ht/DBH.23)
xyplot(final.ht~DBH.23|SPP,data=picea)
xyplot(final.ht~DBH.23|CODE,data=picea)

#-------------------------------------------------------------------------------
#Site Index
#-------------------------------------------------------------------------------
#vicary.site
picea <- picea%>%
  mutate(vicary.si=mapply(vicary.site,SPP="RS",ht=wykoff.ht, age=28))

mean(picea$vicary.si)

picea %>%
  group_by(uid) %>%
  summarise_at(vars(vicary.si), list(name = mean))


#steinman.site
picea <- picea%>%
  mutate(steinman.si=mapply(steinman.site,ht=wykoff.ht, age=28))

mean(picea$steinman.si)

picea %>%
  group_by(uid) %>%
  summarise_at(vars(steinman.si), list(name = mean))

#-------------------------------------------------------------------------------
#Top Height (vicary height)
#-------------------------------------------------------------------------------
picea <- picea%>%
  mutate(topht=mapply(vicary.height,SPP="RS", age=28, si=si))

mean(picea$topht)

picea %>%
  group_by(uid) %>%
  summarise_at(vars(topht), list(name = mean))

#-------------------------------------------------------------------------------
#qmd
#-------------------------------------------------------------------------------
picea <- picea %>%
  mutate(qmd=mapply(qmd,ba=bapa, tpa=tpa))

mean(picea$qmd)

picea %>%
  group_by(uid) %>%
  summarise_at(vars(qmd), list(name = mean))


#-------------------------------------------------------------------------------
#relative density index
#-------------------------------------------------------------------------------
picea <- picea %>%
  mutate(rdi=mapply(relative.density.index,bapa=bapa, qmd=qmd))

mean(picea$rdi)

picea %>%
  group_by(uid) %>%
  summarise_at(vars(rdi), list(name = mean))

xyplot(qmd~bapa|SPP,data=picea)

#-------------------------------------------------------------------------------
#relative spacing
#-------------------------------------------------------------------------------
picea <- picea %>%
  mutate(rs=mapply(relative.spacing,tpa=tpa, domht=topht))

mean(picea$rdi)

picea %>%
  group_by(uid) %>%
  summarise_at(vars(rdi), list(name = mean))

xyplot(qmd~bapa|SPP,data=picea)

#-------------------------------------------------------------------------------
#stand density index
#-------------------------------------------------------------------------------
picea <- picea %>%
  mutate(sdi=mapply(stand.density.index,tpa=tpa, qmd=qmd))

mean(picea$sdi)

picea %>%
  group_by(uid) %>%
  summarise_at(vars(sdi), list(name = mean))

xyplot(qmd~tpa|SPP,data=picea)

#-------------------------------------------------------------------------------
#ANOVA: does diameter vary by species (one-way)
#-------------------------------------------------------------------------------
aov.1<- aov(DBH.23~SPP, data=picea)
summary(aov.1)

TukeyHSD(aov.1)

#Homogeneity of variances
plot(aov.1, 1)

#Normality
plot(aov.1, 2)

#-------------------------------------------------------------------------------
#ANOVA: does diameter vary by species grouped by id, plot level (two-way)
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
#ANOVA: does diameter vary by pure vs mixed (one-way)
#-------------------------------------------------------------------------------
aov.3<- aov(DBH.23~CODE, data=picea)
summary(aov.3)

TukeyHSD(aov.3)

#Homogeneity of variances
plot(aov.3, 1)

#Normality
plot(aov.3, 2)

#-------------------------------------------------------------------------------
#Matt Russell ANOVA Textbook Exercises
#-------------------------------------------------------------------------------
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
            sd.diameter = sd(DBH.23)) %>% 
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

lmodel <- lm(HT.23 ~ DBH.23+SPP, data = picea)
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
#predicting HTs using model need to use heights from intital dataset
#also you have DBH.23 as a mixed and random effect, check Matts code
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

#-------------------------------------------------------------------------------
#ANOVA: does volume vary by species (one-way)
#-------------------------------------------------------------------------------

#Anova
vol.aov.1<- aov(vol~SPP, data=picea)
summary(vol.aov.1)

TukeyHSD(vol.aov.1)

#Homogeneity of variances
plot(vol.aov.1, 1)

#Normality
plot(vol.aov.1, 2)

#or can do similar calculations following Matt Russel's Method
p.vol <- ggplot(picea, aes(factor(SPP), vol)) + 
  geom_boxplot()+
  ylab("Volume (cubic feet)") +
  xlab("Species")
p.vol

vol.aov.2 <- lm(vol ~ SPP, data = picea)

anova(vol.aov.2)

pairwise.t.test(picea$vol, picea$SPP, p.adj = "bonferroni")

library(agricolae)
lsd.picea <- LSD.test(vol.aov.2, "SPP", p.adj = "bonferroni")
lsd.picea$groups

vol.summary <- picea %>% 
  group_by(SPP)  %>%  
  summarize(n.vol = n(),
            mean.vol = mean(vol),
            sd.vol = sd(vol))
vol.summary

limits <- aes(ymax = mean.vol + sd.vol, 
              ymin = mean.vol - sd.vol)
p.vol <- ggplot(vol.summary, aes(SPP, mean.vol)) +
  geom_bar(stat = "identity") +
  geom_errorbar(limits, width = 0.25) +
  ylab("Volume (Cubic Feet") +
  xlab("Species") +
  geom_text(aes(label = c("a", b","c","c","c","c")), vjust = -6) +
  scale_y_continuous(limits = c(0, 5))
p.vol

#-------------------------------------------------------------------------------
#ANOVA: does vol vary by pure vs mixed (one-way)
#-------------------------------------------------------------------------------
p.vol.2 <- ggplot(picea, aes(factor(CODE), vol)) + 
  geom_boxplot()+
  ylab("Volume (cubic feet)") +
  xlab("Species Mix")
p.vol.2
piceamix2.aov <- lm(vol ~ CODE, data = picea)

anova(piceamix2.aov)

pairwise.t.test(picea$vol, picea$CODE, p.adj = "bonferroni")

library(agricolae)
lsd.picea <- LSD.test(piceamix2.aov, "CODE", p.adj = "bonferroni")
lsd.picea$groups

vol.summary2 <- picea %>% 
  group_by(CODE)  %>%  
  summarize(n.vol = n(),
            mean.vol = mean(vol),
            sd.vol = sd(vol))
vol.summary2

limits <- aes(ymax = mean.vol + sd.vol, 
              ymin = mean.vol - sd.vol)
p.vol2 <- ggplot(vol.summary2, aes(CODE, mean.vol)) +
  geom_bar(stat = "identity") +
  geom_errorbar(limits, width = 0.25) +
  ylab("Volume (Cubic Feet)") +
  xlab("Species Mix") +
  geom_text(aes(label = c("a", "ab", "ab", "ab", "abc", "abc", "bcd", "cde", "de", "e")), vjust = -10) +
  scale_y_continuous(limits = c(0, 5))
p.vol2




