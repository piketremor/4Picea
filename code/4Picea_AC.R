#install.packages("devtools")
#install.packages("G:/My Drive/MEForLab",
                 #repos = NULL,
                 #type = "source")
#devtools::install_github("piketremor/MEForLab")

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
picea <- read.csv("4Picea.csv")
stemform <- read.csv("StemForm.csv")

#-------------------------------------------------------------------------------
#join stemform to picea by BLOCK, PLOT, TREE #
#-------------------------------------------------------------------------------
picea$TREE <- as.character(picea$TREE)
stemform$TREE <- as.character(stemform$TREE)
picea <- left_join(picea, stemform, by = c("BLOCK", "PLOT", "TREE", "SPP"))

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
#summary statistics for silviculture lab
#-------------------------------------------------------------------------------
picea$tree.factor <- 10

picea <- picea%>%
  filter(.,STATUS.23=="1")%>%
  mutate(tree.ba = DBH.23^2*0.005454)%>%
  group_by(uid, CODE)%>%
  mutate(bapa = sum(tree.ba*tree.factor), # tree.factor is your expansion factor
  tpa = sum(tree.factor))

plot(picea$bapa,picea$tpa) # most plots are between 20 and 140 bapa

ggplot(picea, aes(x = bapa, y = tpa, color = CODE)) +
  geom_point(alpha = 0.7) +  
  labs(title = "4Picea",
       x = "BAPA",
       y = "TPA",
       color = "CODE") +
  xlim(20, 140) +  
  theme_minimal()

require(lattice)
xyplot(bapa~tpa|CODE,data=picea)

picea$dbh.class <- 2*as.integer((picea$DBH.23+(2/2))/2)


#DBH distribution
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

ggplot(species_proportions, aes(x = interaction(PLOT, BLOCK), y = Proportion, fill = SPP)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "BLOCK and PLOT", y = "Proportion of Trees", title = "Proportion of Each Species by BLOCK and PLOT") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
#-------------------------------------------------------------------------------
#generating height model
#-------------------------------------------------------------------------------
require(nlme)

xyplot(HT.23~DBH.23|SPP,data=picea)
picea$SPP <- as.factor(picea$SPP)
picea$CODE <- as.factor(picea$CODE)
spruce  <- dplyr::filter(picea,SPP=="NS"|SPP=="RS"|SPP=="WS"|SPP=="BS")
xyplot(HT.23~DBH.23|SPP,data=spruce)
str(spruce)


ht.mod <-  nlme(HT.23 ~ 4.5+exp(a+b/(DBH.23+1)),
                data = spruce,
                fixed = a + b ~ 1,
<<<<<<< HEAD
                random = a + b ~ 1 | SPP,  # Random intercept and slope for both
=======
                random = a + b ~ 1 | BLOCK/PLOT/SPP,  
>>>>>>> 69f4450331db5d8e6259a9c46cf6d4be94fbb8c9
                na.action = na.pass,
                start = c(a = 4.5, b = -6),
                control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))
performance(ht.mod)
summary(ht.mod)
ranef(ht.mod)

<<<<<<< HEAD
ht.mod2 <- nlme(HT.23 ~ 4.5+exp(a+b/(DBH.23+1)),
=======
ht.mod <- nlme(HT.23 ~ 4.5 + exp((a + b) / (DBH.23 + 1)),
>>>>>>> 69f4450331db5d8e6259a9c46cf6d4be94fbb8c9
               data = spruce,
               fixed = a + b ~ 1,
               random = a + b ~ 1 | BLOCK/PLOT,  # Random intercept and slope for both
               na.action = na.pass,
               start = c(a = 4.5, b = -6),
               control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))
performance(ht.mod)

AIC(ht.mod,ht.mod2)


ht.mod3 <- nlme(HT.23 ~ 4.5+exp(a+b/(DBH.23+1)),
                data = spruce,
                fixed = a + b ~ SPP,
                random = a + b ~ 1 | BLOCK,  # Random intercept and slope for both
                na.action = na.pass,
                start = c(a = 4.5,0,0,0, b = -6,0,0,0),
                control = nlmeControl(returnObject = TRUE, msMaxIter = 1000000, maxIter = 500000))
                #control=lmeControl(opt="optim"))

AIC(ht.mod3,ht.mod2,ht.mod)
summary(ht.mod3)

<<<<<<< HEAD
ht.mod4 <- nlme(HT.23 ~ 4.5+exp((a+b/DBH.23+1)),
                data = spruce,
                fixed = a + b ~ SPP,
                random = a + b ~ 1 | PLOT,  # Random intercept and slope for both
                na.action = na.pass,
                start = c(a = 4.5,0,0,0, b = -6,0,0,0),
                control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))
summary(ht.mod4)
plot(ht.mod4)
AIC(ht.mod3,ht.mod4,ht.mod2,ht.mod)
=======
ht.mod <- nlme(HT.23 ~ a + b * DBH.23 * Proportion,
                           data = picea,
                           fixed = a + b ~ 1,
                           random = a + b ~ 1 | PLOT/SPP,
                           na.action = na.pass,
                           start = c(a = 4.5, b = -6),
                           control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))
>>>>>>> 69f4450331db5d8e6259a9c46cf6d4be94fbb8c9

ht.mod5 <- nlme(HT.23 ~ 4.5+exp((a+b/DBH.23+1)),
                data = spruce,
                fixed = a + b ~ CODE,
                random = a + b ~ 1 | BLOCK/PLOT,  # Random intercept and slope for both
                na.action = na.pass,
                start = c(a = 4.5,0,0,0,0,0,0,0,0,0,0, b = -6,0,0,0,0,0,0,0,0,0,0),
                control = nlmeControl(returnObject = TRUE, msMaxIter = 10000, maxIter = 5000))

summary(ht.mod5)
AIC(ht.mod,ht.mod2,ht.mod3,ht.mod4,ht.mod5)
spruce$ht.fit <- predict(ht.mod2,spruce)
xyplot(ht.fit~DBH.23|CODE,data=spruce)

d <- (ht.mod2$coefficients$random)
ds <- d$PLOT
BLOCK <- c("B2","B2","B2","B2","B2",
              "B2","B2","B2","B2","B2","B2",
              "B3","B3","B3","B3","B3","B3",
              "B3","B3","B3","B3",
              "B4","B4","B4","B4","B4",
              "B4","B4","B4","B4","B4")
PLOT <- c("1","2","3","4","5","6","7","8","9","10","11",
          "1","2","3","4","5","6","7","8","9","10",
          "1","2","3","4","5","6","7","8","9","10")
re.df <- as.data.frame(ds[,1:2])
red <- cbind(re.df,BLOCK,PLOT)
red$uid <- paste0(red$BLOCK,".",red$PLOT)
red <- red[c(1,2,5)]
spruces <- left_join(spruce,red,by=c("uid"))
pairs(spruces[c(28,29,12,23:24,26)])
plot(spruces$Proportion,spruces$b)
prop.mod <- lm(a~Proportion,data=spruces)
summary(prop.mod)

picea$sp.dum <- ifelse(picea$SPP=="WS"|picea$SPP=="BS"|picea$SPP=="NS"|picea$SPP=="RS",1,0)

#-------------------------------------------------------------------------------
#imputing tree heights (Wykoff)
#-------------------------------------------------------------------------------
require(MEForLab)
# if spruce, use the local model, otherwise use FVS base equations
picea$local.ht <- ifelse(picea$sp.dum=="1",predict(ht.mod2,picea),0)
picea$hd <- picea$local.ht/picea$DBH.23
xyplot(hd~DBH.23|SPP,data=picea,ylim=c(0,20))
picea$local.ht <- ifelse(picea$hd>12,0,picea$local.ht)
xyplot(local.ht~DBH.23|SPP,data=picea)
picea$fvs.ht <- ifelse(picea$local.ht<1,mapply(wykoff.ht,
                       SPP=picea$SPP,
                       DBH=picea$DBH.23),
                       0)
xyplot(fvs.ht~DBH.23|SPP,data=picea)
picea$HT.23[is.na(picea$HT.23)] <- 0
picea$model.ht <- ifelse(picea$local.ht>0,picea$local.ht,picea$fvs.ht)
picea$final.ht <- ifelse(picea$HT.23>0,picea$HT.23,picea$model.ht)


xyplot(final.ht~DBH.23|SPP,data=picea)
xyplot(final.ht~DBH.23|CODE,data=picea)

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

# Filter for the species of interest
picea.ns <- picea %>%
  filter(SPP == "NS")

# Calculate the volume for just NS across plots
picea.ns <- picea.ns %>%
  mutate(vol = mapply(vol.calc, SPP = SPP, DBH = DBH.23, HT = final.ht))

# boxplot
picea.ns<-filter(picea.ns,CODE != "BW"&CODE !="RW"&CODE !="W")
p.vol <- ggplot(picea.ns, aes(factor(CODE), vol)) +
  geom_boxplot()+
  ylab("Volume (cubic feet)") +
  xlab("Species")

p.vol

xyplot(vol~DBH.23|CODE,data=picea)
xyplot(vol~DBH.23|SPP,data=picea)
#-------------------------------------------------------------------------------
#overyielding & transgressive overyielding, looking at the plot level
#-------------------------------------------------------------------------------
total_volume_by_code <- picea %>%
  group_by(BLOCK, PLOT, CODE) %>%
  summarise(Total_Volume = sum(vol, na.rm = TRUE) / 3, .groups = 'drop')

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
<<<<<<< HEAD
  geom_boxplot() +
  labs(title = "Volume by Species Mixture",
       x = "Species Mixture",
       y = "Volume (ft3)") +
  theme_minimal()+
  geom_hline(yintercept=1)
  

=======
  geom_bar(stat = "identity", fill = "grey") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  
  labs(title = "Transgressive Overyielding by Species Mixture", x = "Species Mixture", y = "Transgressive Overyielding Value") +
  theme_minimal()
>>>>>>> 69f4450331db5d8e6259a9c46cf6d4be94fbb8c9

#-------------------------------------------------------------------------------
#Stem Form: volume deductions 
#Remember this is only for height trees, must apply to rest of imputed heights
#-------------------------------------------------------------------------------
# Step 1: generate LR equation to determine what attributes to a defect being present 

picea <- picea %>%
  mutate(defect_present = ifelse( T_DEDUCT > 0, 1, 0))

plot(picea$defect_present)

#NAs make sense since, filter for trees with HTs and Stem Form taken
ggplot(picea, aes(x = as.factor(defect_present))) +
  geom_bar() +
  labs(x = "Defect Present (0 = No, 1 = Yes)", y = "Count") +
  ggtitle("Bar Chart of Defect Presence") +
  theme_minimal()


# Step 2: Defect 1,0 if 1 then go to glm model and account for how much volume is deducted based on SPP, Ht.23, DBH,23, etc. 

#histogram of T_DEDUCT distrbution from 0-100
picea <- picea %>%
  mutate(T_DEDUCT = as.numeric(T_DEDUCT),  # Convert T_DEDUCT to numeric
         defect_present = ifelse(T_DEDUCT > 0, 1, 0))

#keeping only subset of trees with values in T_DEDUCT column 
deduct <- picea %>%
  filter(T_DEDUCT >= 0 & T_DEDUCT <= 100)

ggplot(deduct, aes(x = T_DEDUCT)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(x = "T_DEDUCT", y = "Count") +
  ggtitle("Histogram of T_DEDUCT Values (0-100)") +
  theme_minimal() +
  xlim(0, 100) +  
  ylim(0, 30)  


mod1 <- glm(defect_present ~ DBH.23 + Proportion + HT.23 + HCB.23 + SPP, 
            data = picea, 
            family = binomial(link = "logit"))
summary(mod1)

mod2 <- glm(defect_present ~ DBH.23 + HT.23 + SPP, 
            data = picea, 
            family = binomial(link = "logit"))
summary(mod2)

mod3 <- glm(defect_present ~ DBH.23 + SPP, 
            data = picea, 
            family = binomial(link = "logit"))
summary(mod3)

require(lme4)
mod4 <- glmer(defect_present ~ DBH.23 + (1 |SPP), 
              data = picea, 
              family = binomial(link = "logit"))


AIC(mod1, mod2, mod3, mod4)

# Step 3: use model to to impute T_DEDUCT for all stems

# Step 4: vol - T_DEDUCT = fin.vol

################################ OR ############################################

# Step 1: model with T_DEDUCT on height trees
# Deduct T_DEDUCT from vol where T_DEDUCT is available
picea <- picea %>%
  mutate(vol_adj = ifelse(!is.na(T_DEDUCT), vol * (1 - T_DEDUCT / 100), vol))

# Fit a model to predict T_DEDUCT
mod10 <- glm(T_DEDUCT ~ DBH.23 + SPP, 
                      data = picea, 
                      na.action = na.exclude)  
summary(mod10)

mod11 <- glm(T_DEDUCT ~ DBH.23 * SPP, 
             data = picea, 
             na.action = na.exclude)

mod12 <- glm(T_DEDUCT ~ DBH.23 + HT.23 + SPP + Proportion, 
             data = picea, 
             na.action = na.exclude)  
summary(mod12)

AIC(mod10, mod11, mod12)

# Step 2: apply that model to the rest of the trees
# Predict T_DEDUCT for rows without it
picea$predicted_T_DEDUCT <- predict(mod12, newdata = picea, type = "response")

# Deduct predicted T_DEDUCT from vol for those rows
picea <- picea %>%
  mutate(vol_adj = ifelse(is.na(T_DEDUCT) & !is.na(predicted_T_DEDUCT), 
                               vol * (1 - predicted_T_DEDUCT / 100), 
                               vol_adj))

view(picea)
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
picea$htl[is.na(picea$htl)] <- 0
xyplot(htl~final.ht|SPP,data=picea)
xyplot(htl~final.ht|CODE,data=picea)

#-------------------------------------------------------------------------------
#height/diameter ratios
#-------------------------------------------------------------------------------
# h/d ratio
picea <- picea %>%
  mutate(ht.dbh = final.ht/DBH.23)

xyplot(ht.dbh~DBH.23|SPP,data=picea)
xyplot(ht.dbh~DBH.23|CODE,data=picea)


# omit NAs
picea <- na.omit(picea)


# Fit a mixed-effects model with nested factors
mod1 <- lmer(ht.dbh~log(sp.tpa+0.1)+SPP+(1|PLOT/BLOCK), data=picea)
mod2 <- lmer(ht.dbh~log(sp.bapa+0.1)+SPP+(1|PLOT/BLOCK), data=picea)

# Summary of the model
summary(mod1)
summary(mod2)

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

