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
  labs(x = "PLOT (Grouped by BLOCK)", y = "Proportion of Trees", 
       title = "Proportion of Each Species by BLOCK and PLOT") +
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


plot(picea$bapa, picea$tpa)  # most plots are between 20 and 140 bapa

ggplot(picea, aes(x = bapa, y = tpa, color = CODE)) +
  geom_point(alpha = 0.7) +  
  labs(title = "Picea BAPA vs TPA",
       x = "BAPA",
       y = "TPA",
       color = "CODE") +
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
  labs(title = "Diameter Class Distribution by CODE",
       x = "Diameter Classes (DBH in inches)",
       y = "Number of Trees") +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_blank())  # 0-10" dbh, 0-2, 2-4, 4-6, 6-8, 8-10

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

#ht.mod2 <- nlme(HT.23 ~ 4.5+exp(a+b/(DBH.23+1)),
                #data = picea,
                #fixed = a + b ~ 1,
                #random = a + b ~ 1 | BLOCK/PLOT,  # Random intercept and slope for both
                #na.action = na.pass,
                #start = c(a = 4.5, b = -6),
                #control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))

#ht.mod3 <- nlme(HT.23 ~ 4.5+exp(a+b/(DBH.23+1)),
                #data = spruce,
                #fixed = a + b ~ SPP,
                #random = a + b ~ 1 | BLOCK/PLOT,  # Random intercept and slope for both
                #na.action = na.pass,
                #start = c(a = 4.5,0,0,0, b = -6,0,0,0),
                #control = nlmeControl(returnObject = TRUE, msMaxIter = 1000000, maxIter = 500000))
#control=lmeControl(opt="optim")


#ht.mod4 <- nlme(HT.23 ~ 4.5+exp((a+b/DBH.23+1)),
                #data = spruce,
                #fixed = a + b ~ SPP,
                #random = a + b ~ 1 | PLOT,  # Random intercept and slope for both
                #na.action = na.pass,
                #start = c(a = 4.5,0,0,0, b = -6,0,0,0),
                #control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))

#ht.mod5 <- nlme(HT.23 ~ 4.5+exp((a+b/DBH.23+1)),
                #data = spruce,
                #fixed = a + b ~ CODE,
                #random = a + b ~ 1 | BLOCK/PLOT,  # Random intercept and slope for both
                #na.action = na.pass,
                #start = c(a = 4.5,0,0,0,0,0,0,0,0,0,0, b = -6,0,0,0,0,0,0,0,0,0,0),
                #control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))

#ht.mod6 <- nlme(HT.23 ~ a + b * DBH.23 * Proportion,
                #data = picea,
                #fixed = a + b ~ 1,
                #random = a + b ~ 1 | PLOT/SPP,
                #na.action = na.pass,
                #start = c(a = 4.5, b = -6),
                #control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))
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
#volume calculation at the individual tree level
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

picea$deductclass <- 25*as.integer((picea$T_DEDUCT+(25/2))/25)

response <- "deductclass"  
predictors <- c("SPP", "DBH.23", "HT.23", "HCB.23", "LCR.23", "LAI", "CODE", "elevation", 
                "tri", "tpi", "roughness", "slope", "aspect", "flowdir", "tmin", "tmean", 
                "tmax", "dew", "vpdmax", "vpdmin", "McNab", "Bolstad", "Profile", "Planform", 
                "Winds10", "Winds50", "SWI", "RAD", "ppt", "Parent", "Lithic", "Redox", 
                "dep", "ex.k", "nit", "SWC2", "MeanWD", "Proportion", "Min_depth", "WD2000", 
                "WD2020", "WHC", "ex.mg", "ex.ca", "ph", "bapa", "tpa", "hd", "bal", "htl", 
                "vicary.si", "steinman.si", "topht", "qmd", "rdi", "rs", "sdi", "vol")

picea_varsel <- picea %>%
  filter(!is.na(deductclass) & is.numeric(deductclass)) %>%
  select(c(response, predictors))
picea_varsel[is.na(picea_varsel)] <- 0

# split data into train and test sets
set.seed(123)  # For reproducibility
split <- sample.split(picea_varsel$deductclass, SplitRatio = 0.6666666666666) 
train <- subset(picea_varsel, split == TRUE)
test <- subset(picea_varsel, split == FALSE)

train[] <- lapply(train, function(x) { 
  if(is.numeric(x)) {
    x[is.nan(x)] <- 0
    x[is.infinite(x)] <- 0
  }
  return(x)
})

vsurf <- VSURF(train[, predictors], train[, response], parallel = TRUE, ncores = 3)

vsurf$varselect.pred

names(train[,1:58])

selected_var_cont <- c(bal+sdi+rdi+Min_depth+steinman.si,HT.23)

#-------------------------------------------------------------------------------
#regsubsets variable selection
#-------------------------------------------------------------------------------
library(leaps)

picea_varsel2 <- picea %>%
  filter(!is.na(deductclass) & is.numeric(deductclass))
picea_varsel2[is.na(picea_varsel2)] <- 0


models.iv <- regsubsets(deductclass~SPP + DBH.23 + HT.23 + HCB.23 + LCR.23 + LAI + CODE + elevation + 
                        tri + tpi + roughness + slope + aspect + flowdir + tmin + tmean + 
                        tmax + dew + vpdmax + vpdmin + McNab + Bolstad + Profile + Planform + 
                        Winds10 + Winds50 + SWI + RAD + ppt + Parent + Lithic + Redox + 
                        dep + ex.k + nit + SWC2 + MeanWD + Proportion + Min_depth + WD2000 + 
                        WD2020 + WHC + ex.mg + ex.ca + ph + bapa + tpa + hd + bal+ htl + 
                        vicary.si + steinman.si + topht + qmd + rdi + rs + sdi + vol,data=picea_varsel2,really.big = TRUE,method="exhaustive")

summary(models.iv)
#-------------------------------------------------------------------------------
#volume deduction
#-------------------------------------------------------------------------------
picea2 <- filter(picea,HT.23!="NA")
plot(density(na.omit(picea2$T_DEDUCT))) #look at a zero-inflated model
picea2$T_DEDUCT[is.na(picea2$T_DEDUCT)] <- 0

hist(picea2$T_DEDUCT)
xyplot(T_DEDUCT~DBH.23|CODE,data=picea2) #most deductions in plots with NS, weevil

unique(picea2$T_DEDUCT)
#picea2$deductclass <- 25*as.integer((picea2$T_DEDUCT+(25/2))/25)
hist(picea2$deductclass)

mod1 <- glmer(deductclass~SPP + Min_depth + bal+(1|BLOCK), data=picea2,family="poisson") #mixed poisson
summary(mod1)

mod2 <- glmmTMB(deductclass~SPP+Min_depth+bal+(1|BLOCK), #mixed zero-inflated poisson
                data=picea2,
                ziformula = ~.,
                family=poisson,
                na.action = "na.omit")
summary(mod2)


mod3 <- glmmTMB(deductclass~SPP+Min_depth+bal+(1|BLOCK),
                data=picea2,
                ziformula = ~SPP,
                family=nbinom1,
                na.action = "na.omit")

AIC(mod1,mod2,mod3)


# Matt Russel zero-inflated models
library(pscl)

mod4 <- zeroinfl(deductclass ~ bal+sdi+rdi+Min_depth+steinman.si,
                  data = picea2)
summary(mod4)

mod5 <- zeroinfl(deductclass ~ bal+sdi+rdi+Min_depth+steinman.si, 
                   dist = "negbin", data = picea2)
summary(mod5)



#-------------------------------------------------------------------------------
#overyielding & transgressive overyielding, looking at the plot level
#-------------------------------------------------------------------------------
total_volume_by_code <- picea %>%
  mutate(qmd.2 = qmd(bapa,tpa),
         rd.2 = bapa/sqrt(qmd.2))%>%
  group_by(BLOCK, PLOT, CODE) %>%
  summarise(Total_Volume = sum(vol, na.rm = TRUE) / 3, .groups = 'drop',
            BAPA = mean(bapa),
            TPA = mean(tpa),
            QMD = mean(qmd.2),
            RD = mean(rd.2))

write.csv(total_volume_by_code,"~/Desktop/4Picea_Standlister.csv")


#boxplot of vol by CODE
boxplot_metrics <- total_volume_by_code %>%
  group_by(CODE) %>%
  summarise(
    Q1 = quantile(Total_Volume, 0.25),
    Median = median(Total_Volume),
    Q3 = quantile(Total_Volume, 0.75),
    Min = min(Total_Volume),
    Max = max(Total_Volume),
    .groups = 'drop'
  )


ggplot(total_volume_by_code, aes(x = CODE, y = Total_Volume)) +
  geom_boxplot() +
  labs(title = "Volume by Species Mixture",
       x = "Species Mixture",
       y = "Volume (ft3)") +
  theme_minimal()

data <- read.table(text = "
CODE    VOL
NW	    57
BW	    49
RW	    40
W	      51
NR	    48
BN	    45
N	      56
BR	    35
R	      17
B	      42
 ", header = TRUE, stringsAsFactors = FALSE)


# Function to calculate overyielding

calculate_overyielding <- function(data) {
  results <- data.frame(Mixture = character(), Overyielding = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:nrow(data)) {
    mixture <- data$CODE[i]
    
    # Check if the mixture consists of two species (like 'NW')
    if (nchar(mixture) == 2) {
      species1 <- substr(mixture, 1, 1)  # First species
      species2 <- substr(mixture, 2, 2)  # Second species
      
      # Get yields for the two monocultures
      volume1 <- data$VOL[data$CODE == species1]
      volume2 <- data$VOL[data$CODE == species2]
      
      # Get the mixture yield
      mixture_volume <- data$VOL[i]
      
      # Calculate the average monoculture volume
      if (length(volume1) > 0 && length(volume2) > 0) {
        average_volume <- mean(c(volume1, volume2), na.rm = TRUE)
        overyielding <- mixture_volume / average_volume
        
        # Store the results
        results <- rbind(results, data.frame(Mixture = mixture, Overyielding = overyielding))
      }
    }
  }
  
  return(results)
}

overyielding_results <- calculate_overyielding(data)
print(overyielding_results)

ggplot(overyielding_results, aes(x = Mixture, y = Overyielding)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Reference line at overyielding = 1
  labs(title = "Overyielding by Species Mixture", x = "Species Mixture", y = "Overyielding Value") +
  theme_minimal()

# Function to calculate transgressive overyielding
calculate_transgressive_overyielding <- function(data) {
  results <- data.frame(Mixture = character(), Transgressive_Overyielding = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:nrow(data)) {
    mixture <- data$CODE[i]
    
    # Check if the mixture consists of two species (like 'NW')
    if (nchar(mixture) == 2) {
      species1 <- substr(mixture, 1, 1)  # First species
      species2 <- substr(mixture, 2, 2)  # Second species
      
      # Get yields for the two monocultures
      volume1 <- data$VOL[data$CODE == species1]
      volume2 <- data$VOL[data$CODE == species2]
      
      # Get the mixture yield
      mixture_volume <- data$VOL[i]
      
      # Calculate the maximum monoculture volume
      if (length(volume1) > 0 && length(volume2) > 0) {
        max_volume <- max(volume1, volume2, na.rm = TRUE)
        transgressive_overyielding <- mixture_volume / max_volume
        
        # Store the results
        results <- rbind(results, data.frame(Mixture = mixture, Transgressive_Overyielding = transgressive_overyielding))
      }
    }
  }
  
  return(results)
}

transgressive_overyielding_results <- calculate_transgressive_overyielding(data)
print(transgressive_overyielding_results)

ggplot(transgressive_overyielding_results, aes(x = Mixture, y = Transgressive_Overyielding)) +
  geom_bar(stat = "identity", fill = "grey") +  # Use geom_bar for bar chart
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Reference line at 1
  labs(title = "Transgressive Overyielding by Species Mixture", 
       x = "Species Mixture", 
       y = "Transgressive Overyielding Value") +
  theme_minimal()
