# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 11 - Calibration plots
# Author: Nathan D. Hooven
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Email: nathan.hooven@uky.edu
# Date started: 12 May 2022
# Date completed: 13 May 2022
# Date modified: 
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(survival)      # conditional logistic regressions
library(mefa4)         # %notin% function
library(mvtnorm)       # multivariate normal distribution
library(KernSmooth)    # calculate simulation envelopes
library(cowplot)       # arranging multiple plots

#__________________________________________________________________________________________________
# 2. Read in data ----
#__________________________________________________________________________________________________

load("SSF_models.RData")

#_____________________________________________________________________________________________________________
# 3. Split strata into test and training sets (3 folds) ----
#_____________________________________________________________________________________________________________
# 3a. Spring ----
#_____________________________________________________________________________________________________________

# determine number of strata, divide by two and round up
n.test.spring <- round(length(unique(deer.data.spring.s$step_id_)) / 2)

# sample strata used locations for testing set
test.spring <- deer.data.spring.s %>% filter(case_ == 1) %>% 
                                      slice_sample(n = n.test.spring) 

# bring random locations along for the full test set
test.spring.all <- deer.data.spring.s %>% filter(step_id_ %in% unique(test.spring$step_id_))

# use the other strata as the training set
train.spring.all <- deer.data.spring.s %>% filter(step_id_ %notin% unique(test.spring$step_id_))

#_____________________________________________________________________________________________________________
# 3b. Summer ----
#_____________________________________________________________________________________________________________

# determine number of strata, divide by two and round up
n.test.summer <- round(length(unique(deer.data.summer.s$step_id_)) / 2)

# sample strata used locations for testing set
test.summer <- deer.data.summer.s %>% filter(case_ == 1) %>% 
                                      slice_sample(n = n.test.summer) 

# bring random locations along for the full test set
test.summer.all <- deer.data.summer.s %>% filter(step_id_ %in% unique(test.summer$step_id_))

# use the other strata as the training set
train.summer.all <- deer.data.summer.s %>% filter(step_id_ %notin% unique(test.summer$step_id_))

#_____________________________________________________________________________________________________________
# 3c. Fall ----
#_____________________________________________________________________________________________________________

# determine number of strata, divide by two and round up
n.test.fall <- round(length(unique(deer.data.fall.s$step_id_)) / 2)

# sample strata used locations for testing set
test.fall <- deer.data.fall.s %>% filter(case_ == 1) %>% 
                                      slice_sample(n = n.test.fall) 

# bring random locations along for the full test set
test.fall.all <- deer.data.fall.s %>% filter(step_id_ %in% unique(test.fall$step_id_))

# use the other strata as the training set
train.fall.all <- deer.data.fall.s %>% filter(step_id_ %notin% unique(test.fall$step_id_))

#_____________________________________________________________________________________________________________
# 3d. Dispersal ----
#_____________________________________________________________________________________________________________

# determine number of strata, divide by two and round up
n.test.disp <- round(length(unique(deer.data.disp.s$step_id_)) / 2)

# sample strata used locations for testing set
test.disp <- deer.data.disp.s %>% filter(case_ == 1) %>% 
                                      slice_sample(n = n.test.disp) 

# bring random locations along for the full test set
test.disp.all <- deer.data.disp.s %>% filter(step_id_ %in% unique(test.disp$step_id_))

# use the other strata as the training set
train.disp.all <- deer.data.disp.s %>% filter(step_id_ %notin% unique(test.disp$step_id_))

#_____________________________________________________________________________________________________________
# 4. Fit models to training data ----
#_____________________________________________________________________________________________________________
# 4a. Spring ----
#_____________________________________________________________________________________________________________

spring.train.mod.1 <- clogit(case_ ~
                               dForest.P +
                               dStream + dStream.Q +
                               strata(step_id_),
                             data = train.spring.all)

spring.train.mod.2 <- clogit(case_ ~
                               dForest.P +
                               dStream + dStream.Q +
                               D2.all.hr.175.P +
                               strata(step_id_),
                             data = train.spring.all)

spring.train.mod.3 <- clogit(case_ ~
                               dForest.P +
                               dStream + dStream.Q +
                                D2.all.hr.175.P +
                               D2.indiv.hr.500 + D2.indiv.hr.500.Q +
                               strata(step_id_),
                             data = train.spring.all)

#_____________________________________________________________________________________________________________
# 4b. Summer ----
#_____________________________________________________________________________________________________________

summer.train.mod.1 <- clogit(case_ ~
                               dForest.P +
                               dStream +
                               strata(step_id_),
                             data = train.summer.all)

summer.train.mod.2 <- clogit(case_ ~
                               dForest.P +
                               dStream +
                               D2.all.hr.250.P +
                               strata(step_id_),
                             data = train.summer.all)

summer.train.mod.3 <- clogit(case_ ~
                               dForest.P +
                               dStream +
                               D2.all.hr.250.P +
                               D2.indiv.hr.175.P + 
                               strata(step_id_),
                             data = train.summer.all)

#_____________________________________________________________________________________________________________
# 4c. Fall ----
#_____________________________________________________________________________________________________________

fall.train.mod.1 <- clogit(case_ ~
                               dAg.P +
                                dForest.P +
                                dRoad.P +
                               dStream +
                               strata(step_id_),
                             data = train.fall.all)

fall.train.mod.2 <- clogit(case_ ~
                               dAg.P +
                                dForest.P +
                                dRoad.P +
                               dStream +
                               D2.all.hr.175.P +
                               strata(step_id_),
                             data = train.fall.all)

fall.train.mod.3 <- clogit(case_ ~
                               dAg.P +
                                dForest.P +
                                dRoad.P +
                               dStream +
                               D2.all.hr.175.P +
                               D2.indiv.hr.350.P +
                               strata(step_id_),
                             data = train.fall.all)

#_____________________________________________________________________________________________________________
# 4d. Dispersal ----
#_____________________________________________________________________________________________________________

disp.train.mod.1 <- clogit(case_ ~
                               dForest.P +
                               dStream.P +
                               strata(step_id_),
                             data = train.disp.all)

disp.train.mod.2 <- clogit(case_ ~
                               dForest.P +
                               dStream.P +
                               D2.all.hr.175.P +
                               strata(step_id_),
                             data = train.disp.all)

disp.train.mod.3 <- clogit(case_ ~
                               dForest.P +
                               dStream.P +
                               D2.all.hr.175.P +
                               D2.indiv.hr.175.P +
                               strata(step_id_),
                             data = train.disp.all)

#_____________________________________________________________________________________________________________
# 5. Draw possible parameter values from a multivariate normal distribution ----
#_____________________________________________________________________________________________________________
# 5a. Spring ----
#_____________________________________________________________________________________________________________

# model 1
spring.train.mod.1.param <- spring.train.mod.1$coefficients

spring.train.mod.1.vcov <- vcov(spring.train.mod.1)

spring.train.mod.1.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = spring.train.mod.1.param, 
                                                sigma = spring.train.mod.1.vcov))

# model 2
spring.train.mod.2.param <- spring.train.mod.2$coefficients

spring.train.mod.2.vcov <- vcov(spring.train.mod.2)

spring.train.mod.2.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = spring.train.mod.2.param, 
                                                sigma = spring.train.mod.2.vcov))

# model 3
spring.train.mod.3.param <- spring.train.mod.3$coefficients

spring.train.mod.3.vcov <- vcov(spring.train.mod.3)

spring.train.mod.3.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = spring.train.mod.3.param, 
                                                sigma = spring.train.mod.3.vcov))

#_____________________________________________________________________________________________________________
# 5b. Summer ----
#_____________________________________________________________________________________________________________

# model 1
summer.train.mod.1.param <- summer.train.mod.1$coefficients

summer.train.mod.1.vcov <- vcov(summer.train.mod.1)

summer.train.mod.1.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = summer.train.mod.1.param, 
                                                sigma = summer.train.mod.1.vcov))

# model 2
summer.train.mod.2.param <- summer.train.mod.2$coefficients

summer.train.mod.2.vcov <- vcov(summer.train.mod.2)

summer.train.mod.2.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = summer.train.mod.2.param, 
                                                sigma = summer.train.mod.2.vcov))

# model 3
summer.train.mod.3.param <- summer.train.mod.3$coefficients

summer.train.mod.3.vcov <- vcov(summer.train.mod.3)

summer.train.mod.3.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = summer.train.mod.3.param, 
                                                sigma = summer.train.mod.3.vcov))

#_____________________________________________________________________________________________________________
# 5c. Fall ----
#_____________________________________________________________________________________________________________

# model 1
fall.train.mod.1.param <- fall.train.mod.1$coefficients

fall.train.mod.1.vcov <- vcov(fall.train.mod.1)

fall.train.mod.1.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = fall.train.mod.1.param, 
                                                sigma = fall.train.mod.1.vcov))

# model 2
fall.train.mod.2.param <- fall.train.mod.2$coefficients

fall.train.mod.2.vcov <- vcov(fall.train.mod.2)

fall.train.mod.2.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = fall.train.mod.2.param, 
                                                sigma = fall.train.mod.2.vcov))

# model 3
fall.train.mod.3.param <- fall.train.mod.3$coefficients

fall.train.mod.3.vcov <- vcov(fall.train.mod.3)

fall.train.mod.3.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = fall.train.mod.3.param, 
                                                sigma = fall.train.mod.3.vcov))

#_____________________________________________________________________________________________________________
# 5d. Dispersal ----
#_____________________________________________________________________________________________________________

# model 1
disp.train.mod.1.param <- disp.train.mod.1$coefficients

disp.train.mod.1.vcov <- vcov(disp.train.mod.1)

disp.train.mod.1.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = disp.train.mod.1.param, 
                                                sigma = disp.train.mod.1.vcov))

# model 2
disp.train.mod.2.param <- disp.train.mod.2$coefficients

disp.train.mod.2.vcov <- vcov(disp.train.mod.2)

disp.train.mod.2.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = disp.train.mod.2.param, 
                                                sigma = disp.train.mod.2.vcov))

# model 3
disp.train.mod.3.param <- disp.train.mod.3$coefficients

disp.train.mod.3.vcov <- vcov(disp.train.mod.3)

disp.train.mod.3.mvt <- as.data.frame(rmvnorm(n = 500, 
                                                mean = disp.train.mod.3.param, 
                                                sigma = disp.train.mod.3.vcov))

#_____________________________________________________________________________________________________________
# 6. Calculate RSF scores ----
#_____________________________________________________________________________________________________________
# 6a. Spring ----
#_____________________________________________________________________________________________________________

# model 1
spring.rsf.scores.1 <- data.frame("stratum" = test.spring.all$step_id_)

for (y in 1:500) {
  
  spring.rsf.score.1 <- data.frame(exp(spring.train.mod.1.mvt[y, 1]*test.spring.all$dForest.P +
                                       spring.train.mod.1.mvt[y, 2]*test.spring.all$dStream +
                                       spring.train.mod.1.mvt[y, 3]*test.spring.all$dStream.Q))
  
  spring.rsf.scores.1 <- cbind(spring.rsf.scores.1, spring.rsf.score.1)
  
}

# change column names
colnames(spring.rsf.scores.1) <- c("stratum", 1:500)

# model 2
spring.rsf.scores.2 <- data.frame("stratum" = test.spring.all$step_id_)

for (y in 1:500) {
  
  spring.rsf.score.2 <- data.frame(exp(spring.train.mod.2.mvt[y, 1]*test.spring.all$dForest.P +
                                       spring.train.mod.2.mvt[y, 2]*test.spring.all$dStream +
                                       spring.train.mod.2.mvt[y, 3]*test.spring.all$dStream.Q +
                                       spring.train.mod.2.mvt[y, 4]*test.spring.all$D2.all.hr.175.P))
  
  spring.rsf.scores.2 <- cbind(spring.rsf.scores.2, spring.rsf.score.2)
  
}

# change column names
colnames(spring.rsf.scores.2) <- c("stratum", 1:500)

# model 3
spring.rsf.scores.3 <- data.frame("stratum" = test.spring.all$step_id_)

for (y in 1:500) {
  
  spring.rsf.score.3 <- data.frame(exp(spring.train.mod.3.mvt[y, 1]*test.spring.all$dForest.P +
                                       spring.train.mod.3.mvt[y, 2]*test.spring.all$dStream +
                                       spring.train.mod.3.mvt[y, 3]*test.spring.all$dStream.Q +
                                       spring.train.mod.3.mvt[y, 4]*test.spring.all$D2.all.hr.175.P +
                                       spring.train.mod.3.mvt[y, 5]*test.spring.all$D2.indiv.hr.500 +
                                       spring.train.mod.3.mvt[y, 6]*test.spring.all$D2.indiv.hr.500.Q))
  
  spring.rsf.scores.3 <- cbind(spring.rsf.scores.3, spring.rsf.score.3)
  
}

# change column names
colnames(spring.rsf.scores.3) <- c("stratum", 1:500)

#_____________________________________________________________________________________________________________
# 6b. Summer ----
#_____________________________________________________________________________________________________________

# model 1
summer.rsf.scores.1 <- data.frame("stratum" = test.summer.all$step_id_)

for (y in 1:500) {
  
  summer.rsf.score.1 <- data.frame(exp(summer.train.mod.1.mvt[y, 1]*test.summer.all$dForest.P +
                                       summer.train.mod.1.mvt[y, 2]*test.summer.all$dStream))
  
  summer.rsf.scores.1 <- cbind(summer.rsf.scores.1, summer.rsf.score.1)
  
}

# change column names
colnames(summer.rsf.scores.1) <- c("stratum", 1:500)

# model 2
summer.rsf.scores.2 <- data.frame("stratum" = test.summer.all$step_id_)

for (y in 1:500) {
  
  summer.rsf.score.2 <- data.frame(exp(summer.train.mod.2.mvt[y, 1]*test.summer.all$dForest.P +
                                       summer.train.mod.2.mvt[y, 2]*test.summer.all$dStream +
                                       summer.train.mod.2.mvt[y, 3]*test.summer.all$D2.all.hr.250.P))
  
  summer.rsf.scores.2 <- cbind(summer.rsf.scores.2, summer.rsf.score.2)
  
}

# change column names
colnames(summer.rsf.scores.2) <- c("stratum", 1:500)

# model 3
summer.rsf.scores.3 <- data.frame("stratum" = test.summer.all$step_id_)

for (y in 1:500) {
  
  summer.rsf.score.3 <- data.frame(exp(summer.train.mod.3.mvt[y, 1]*test.summer.all$dForest.P +
                                       summer.train.mod.3.mvt[y, 2]*test.summer.all$dStream +
                                       summer.train.mod.3.mvt[y, 3]*test.summer.all$D2.all.hr.250.P +
                                       summer.train.mod.3.mvt[y, 4]*test.summer.all$D2.indiv.hr.175.P))
  
  summer.rsf.scores.3 <- cbind(summer.rsf.scores.3, summer.rsf.score.3)
  
}

# change column names
colnames(summer.rsf.scores.3) <- c("stratum", 1:500)

#_____________________________________________________________________________________________________________
# 6c. Fall ----
#_____________________________________________________________________________________________________________

# model 1
fall.rsf.scores.1 <- data.frame("stratum" = test.fall.all$step_id_)

for (y in 1:500) {
  
  fall.rsf.score.1 <- data.frame(exp(fall.train.mod.1.mvt[y, 1]*test.fall.all$dAg.P +
                                     fall.train.mod.1.mvt[y, 2]*test.fall.all$dForest.P +
                                     fall.train.mod.1.mvt[y, 3]*test.fall.all$dRoad.P +
                                     fall.train.mod.1.mvt[y, 4]*test.fall.all$dStream))
  
  fall.rsf.scores.1 <- cbind(fall.rsf.scores.1, fall.rsf.score.1)
  
}

# change column names
colnames(fall.rsf.scores.1) <- c("stratum", 1:500)

# model 2
fall.rsf.scores.2 <- data.frame("stratum" = test.fall.all$step_id_)

for (y in 1:500) {
  
  fall.rsf.score.2 <- data.frame(exp(fall.train.mod.2.mvt[y, 1]*test.fall.all$dAg.P +
                                     fall.train.mod.2.mvt[y, 2]*test.fall.all$dForest.P +
                                     fall.train.mod.2.mvt[y, 3]*test.fall.all$dRoad.P +
                                     fall.train.mod.2.mvt[y, 4]*test.fall.all$dStream +
                                     fall.train.mod.2.mvt[y, 5]*test.fall.all$D2.all.hr.175.P))
  
  fall.rsf.scores.2 <- cbind(fall.rsf.scores.2, fall.rsf.score.2)
  
}

# change column names
colnames(fall.rsf.scores.2) <- c("stratum", 1:500)

# model 3
fall.rsf.scores.3 <- data.frame("stratum" = test.fall.all$step_id_)

for (y in 1:500) {
  
  fall.rsf.score.3 <- data.frame(exp(fall.train.mod.3.mvt[y, 1]*test.fall.all$dAg.P +
                                     fall.train.mod.3.mvt[y, 2]*test.fall.all$dForest.P +
                                     fall.train.mod.3.mvt[y, 3]*test.fall.all$dRoad.P +
                                     fall.train.mod.3.mvt[y, 4]*test.fall.all$dStream +
                                     fall.train.mod.3.mvt[y, 5]*test.fall.all$D2.all.hr.175.P +
                                     fall.train.mod.3.mvt[y, 6]*test.fall.all$D2.indiv.hr.350.P))
  
  fall.rsf.scores.3 <- cbind(fall.rsf.scores.3, fall.rsf.score.3)
  
}

# change column names
colnames(fall.rsf.scores.3) <- c("stratum", 1:500)

#_____________________________________________________________________________________________________________
# 6d. Dispersal ----
#_____________________________________________________________________________________________________________

# model 1
disp.rsf.scores.1 <- data.frame("stratum" = test.disp.all$step_id_)

for (y in 1:500) {
  
  disp.rsf.score.1 <- data.frame(exp(disp.train.mod.1.mvt[y, 1]*test.disp.all$dForest.P +
                                     disp.train.mod.1.mvt[y, 2]*test.disp.all$dStream.P))
  
  disp.rsf.scores.1 <- cbind(disp.rsf.scores.1, disp.rsf.score.1)
  
}

# change column names
colnames(disp.rsf.scores.1) <- c("stratum", 1:500)

# model 2
disp.rsf.scores.2 <- data.frame("stratum" = test.disp.all$step_id_)

for (y in 1:500) {
  
  disp.rsf.score.2 <- data.frame(exp(disp.train.mod.2.mvt[y, 1]*test.disp.all$dForest.P +
                                     disp.train.mod.2.mvt[y, 2]*test.disp.all$dStream.P +
                                     disp.train.mod.2.mvt[y, 3]*test.disp.all$D2.all.hr.175.P))
  
  disp.rsf.scores.2 <- cbind(disp.rsf.scores.2, disp.rsf.score.2)
  
}

# change column names
colnames(disp.rsf.scores.2) <- c("stratum", 1:500)

# model 3
disp.rsf.scores.3 <- data.frame("stratum" = test.disp.all$step_id_)

for (y in 1:500) {
  
  disp.rsf.score.3 <- data.frame(exp(disp.train.mod.3.mvt[y, 1]*test.disp.all$dForest.P +
                                     disp.train.mod.3.mvt[y, 2]*test.disp.all$dStream.P +
                                     disp.train.mod.3.mvt[y, 3]*test.disp.all$D2.all.hr.175.P +
                                     disp.train.mod.3.mvt[y, 4]*test.disp.all$D2.indiv.hr.175.P))
  
  disp.rsf.scores.3 <- cbind(disp.rsf.scores.3, disp.rsf.score.3)
  
}

# change column names
colnames(disp.rsf.scores.3) <- c("stratum", 1:500)

#_____________________________________________________________________________________________________________
# 7. Sample from test data (one from each stratum) ----
#_____________________________________________________________________________________________________________
# 7a. Spring ----
#_____________________________________________________________________________________________________________

# model 1
sampled.rows.spring.1 <- data.frame()

for (y in 108:500) {
  
  value <- y
  
  for (x in unique(test.spring.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.spring.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = spring.rsf.scores.1[spring.rsf.scores.1$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.spring.1 <- rbind(sampled.rows.spring.1, sampled.row)
    
  }
  
}

# model 2
sampled.rows.spring.2 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(test.spring.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.spring.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = spring.rsf.scores.2[spring.rsf.scores.2$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.spring.2 <- rbind(sampled.rows.spring.2, sampled.row)
    
  }
  
}

# model 3
sampled.rows.spring.3 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(test.spring.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.spring.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = spring.rsf.scores.3[spring.rsf.scores.3$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.spring.3 <- rbind(sampled.rows.spring.3, sampled.row)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7b. Summer ----
#_____________________________________________________________________________________________________________

# model 1
sampled.rows.summer.1 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(test.summer.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.summer.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = summer.rsf.scores.1[summer.rsf.scores.1$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.summer.1 <- rbind(sampled.rows.summer.1, sampled.row)
    
  }
  
}

# model 2
sampled.rows.summer.2 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(test.summer.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.summer.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = summer.rsf.scores.2[summer.rsf.scores.2$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.summer.2 <- rbind(sampled.rows.summer.2, sampled.row)
    
  }
  
}

# model 3
sampled.rows.summer.3 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(test.summer.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.summer.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = summer.rsf.scores.3[summer.rsf.scores.3$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.summer.3 <- rbind(sampled.rows.summer.3, sampled.row)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7c. Fall ----
#_____________________________________________________________________________________________________________

# model 1
sampled.rows.fall.1 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(test.fall.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.fall.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = fall.rsf.scores.1[fall.rsf.scores.1$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.fall.1 <- rbind(sampled.rows.fall.1, sampled.row)
    
  }
  
}

# model 2
sampled.rows.fall.2 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(test.fall.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.fall.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = fall.rsf.scores.2[fall.rsf.scores.2$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.fall.2 <- rbind(sampled.rows.fall.2, sampled.row)
    
  }
  
}

# model 3
sampled.rows.fall.3 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(test.fall.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.fall.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = fall.rsf.scores.3[fall.rsf.scores.3$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.fall.3 <- rbind(sampled.rows.fall.3, sampled.row)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7d. Dispersal ----
#_____________________________________________________________________________________________________________

# model 1
sampled.rows.disp.1 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(test.disp.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.disp.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = disp.rsf.scores.1[disp.rsf.scores.1$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.disp.1 <- rbind(sampled.rows.disp.1, sampled.row)
    
  }
  
}

# model 2
sampled.rows.disp.2 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(test.disp.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.disp.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = disp.rsf.scores.2[disp.rsf.scores.2$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.disp.2 <- rbind(sampled.rows.disp.2, sampled.row)
    
  }
  
}

# model 3
sampled.rows.disp.3 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(test.disp.all$step_id_)) {
    
    stratum <- x
    
    stratum.data <- test.disp.all %>% filter(step_id_ == stratum) %>%
      bind_cols(RSF.score = disp.rsf.scores.3[disp.rsf.scores.3$stratum == stratum, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
   sampled.rows.disp.3 <- rbind(sampled.rows.disp.3, sampled.row)
    
  }
  
}

save.image("sampled_uhc_5_12_22.RData")

#_____________________________________________________________________________________________________________
# 8. Create simulation envelopes ----
#_____________________________________________________________________________________________________________
# 8a. Write function ----
#_____________________________________________________________________________________________________________

sim_env <- function(variable, season, x.val = 500) {
  
  # create blank df
  sim.env <- data.frame()
  
  # define input datasets depending upon season
  if (season == "spring") {
    
    # range of values
    range.val.1 <- c(min(sampled.rows.spring.1[variable]), max(sampled.rows.spring.1[variable]))
    range.val.2 <- c(min(sampled.rows.spring.2[variable]), max(sampled.rows.spring.2[variable]))
    range.val.3 <- c(min(sampled.rows.spring.3[variable]), max(sampled.rows.spring.3[variable]))
    
    # sampled rows
    sampled.rows.1 <- sampled.rows.spring.1
    sampled.rows.2 <- sampled.rows.spring.2
    sampled.rows.3 <- sampled.rows.spring.3
    
  } else { 
    
    if (season == "summer") {
    
    # range of values
    range.val.1 <- c(min(sampled.rows.summer.1[variable]), max(sampled.rows.summer.1[variable]))
    range.val.2 <- c(min(sampled.rows.summer.2[variable]), max(sampled.rows.summer.2[variable]))
    range.val.3 <- c(min(sampled.rows.summer.3[variable]), max(sampled.rows.summer.3[variable]))
    
    # sampled rows
    sampled.rows.1 <- sampled.rows.summer.1
    sampled.rows.2 <- sampled.rows.summer.2
    sampled.rows.3 <- sampled.rows.summer.3
    
  } else {
    
    if (season == "fall") {
    
    # range of values
    range.val.1 <- c(min(sampled.rows.fall.1[variable]), max(sampled.rows.fall.1[variable]))
    range.val.2 <- c(min(sampled.rows.fall.2[variable]), max(sampled.rows.fall.2[variable]))
    range.val.3 <- c(min(sampled.rows.fall.3[variable]), max(sampled.rows.fall.3[variable]))
    
    # sampled rows
    sampled.rows.1 <- sampled.rows.fall.1
    sampled.rows.2 <- sampled.rows.fall.2
    sampled.rows.3 <- sampled.rows.fall.3
    
  } else { 
    
    if (season == "disp") {
    
    # range of values
    range.val.1 <- c(min(sampled.rows.disp.1[variable]), max(sampled.rows.disp.1[variable]))
    range.val.2 <- c(min(sampled.rows.disp.2[variable]), max(sampled.rows.disp.2[variable]))
    range.val.3 <- c(min(sampled.rows.disp.3[variable]), max(sampled.rows.disp.3[variable]))
    
    # sampled rows
    sampled.rows.1 <- sampled.rows.disp.1
    sampled.rows.2 <- sampled.rows.disp.2
    sampled.rows.3 <- sampled.rows.disp.3 
    
    }
  }
  }
  }
    
  # model 1
  # how many x values?
  dens.1 <- matrix(NA, x.val, x.val)
  
  # for loop to iterate over 500 possible values
  for (z in 1:500) {
    
    which.run <- z
    
    run.rows <- sampled.rows.1 %>% filter(run == which.run)
    
    run.dens <- bkde(c(run.rows[variable][,1]), range.x = range.val.1, gridsize = 500)
    
    dens.1[z, ] <- run.dens$y
    
  }
  
  # compute 95% simulation envelope
  mean.1 <- apply(dens.1, 2, mean, na.rm = TRUE)
  
  low.1 <- apply(dens.1, 2, quantile, prob = 0.025, na.rm = TRUE)
  
  high.1 <- apply(dens.1, 2, quantile, prob = 0.975, na.rm = TRUE)
  
  # create a df
  sim.env.1 <- data.frame("var" = variable,
                          "model" = 1,
                          "mean" = mean.1,
                          "low" = low.1,
                          "high" = high.1,
                          "x" = seq(min(sampled.rows.1[variable]), 
                                    max(sampled.rows.1[variable]), 
                                    length.out = x.val))
  
  # bind to sim.env df
  sim.env <- rbind(sim.env, sim.env.1)
  
  # model 2
  # how many x values?
  dens.2 <- matrix(NA, x.val, x.val)
  
  # for loop to iterate over 500 possible values
  for (z in 1:500) {
    
    which.run <- z
    
    run.rows <- sampled.rows.2 %>% filter(run == which.run)
    
    run.dens <- bkde(c(run.rows[variable][ ,1]), range.x = range.val.2, gridsize = 500)
    
    dens.2[z, ] <- run.dens$y
    
  }
  
  # compute 95% simulation envelope
  mean.2 <- apply(dens.2, 2, mean, na.rm = TRUE)
  
  low.2 <- apply(dens.2, 2, quantile, prob = 0.025, na.rm = TRUE)
  
  high.2 <- apply(dens.2, 2, quantile, prob = 0.975, na.rm = TRUE)
  
  # create a df
  sim.env.2 <- data.frame("var" = variable,
                          "model" = 2,
                          "mean" = mean.2,
                          "low" = low.2,
                          "high" = high.2,
                          "x" = seq(min(sampled.rows.2[variable]), 
                                    max(sampled.rows.2[variable]), 
                                    length.out = x.val))
  
  # bind to sim.env df
  sim.env <- rbind(sim.env, sim.env.2)
  
  # model 3
  # how many x values?
  dens.3 <- matrix(NA, x.val, x.val)
  
  # for loop to iterate over 500 possible values
  for (z in 1:500) {
    
    which.run <- z
    
    run.rows <- sampled.rows.3 %>% filter(run == which.run)
    
    run.dens <- bkde(c(run.rows[variable][ ,1]), range.x = range.val.3, gridsize = 500)
    
    dens.3[z, ] <- run.dens$y
    
  }
  
  # compute 95% simulation envelope
  mean.3 <- apply(dens.3, 2, mean, na.rm = TRUE)
  
  low.3 <- apply(dens.3, 2, quantile, prob = 0.025, na.rm = TRUE)
  
  high.3 <- apply(dens.3, 2, quantile, prob = 0.975, na.rm = TRUE)
  
  # create a df
  sim.env.3 <- data.frame("var" = variable,
                          "model" = 3,
                          "mean" = mean.3,
                          "low" = low.3,
                          "high" = high.3,
                          "x" = seq(min(sampled.rows.3[variable]), 
                                    max(sampled.rows.3[variable]), 
                                    length.out = x.val))
  
  # bind to sim.env df
  sim.env <- rbind(sim.env, sim.env.3)
  
  # return the final df
  return(sim.env)
  
}

#_____________________________________________________________________________________________________________
# 8. Simulation envelopes per variable ----
#_____________________________________________________________________________________________________________

# should we just include the linear terms here? does it make sense to use the log variables?

sim.spring <- sim_env(variable = "D2.indiv.hr.500",
                      season = "spring")

sim.summer <- sim_env(variable = "D2.indiv.hr.175.P",
                      season = "summer")

sim.fall <- sim_env(variable = "D2.indiv.hr.350.P",
                      season = "fall")

sim.disp <- sim_env(variable = "D2.indiv.hr.175.P",
                      season = "disp")

#_____________________________________________________________________________________________________________
# 9. Used-habitat calibration plots ----
#_____________________________________________________________________________________________________________
# 9a. Spring ----
#_____________________________________________________________________________________________________________

uhc.spring <- ggplot(data = sim.spring) +
       theme_bw() +
       facet_wrap(~model) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = test.spring.all[test.spring.all$case_ == 0,], 
                    aes(x = D2.indiv.hr.500),
                    linetype = "dashed") +
       geom_density(data = test.spring.all[test.spring.all$case_ == 1,], 
                    aes(x = D2.indiv.hr.500),
                    size = 1.25,
                    color = "darkgreen") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank(),
             axis.title.x = element_blank()) +
       xlab(expression(scale(D_IND[500]))) +
       coord_cartesian(xlim = c(-0.5, 2))

#_____________________________________________________________________________________________________________
# 9b. Summer ----
#_____________________________________________________________________________________________________________

uhc.summer <- ggplot(data = sim.summer) +
       theme_bw() +
       facet_wrap(~model) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = test.summer.all[test.summer.all$case_ == 0,], 
                    aes(x = D2.indiv.hr.175.P),
                    linetype = "dashed") +
       geom_density(data = test.summer.all[test.summer.all$case_ == 1,], 
                    aes(x = D2.indiv.hr.175.P),
                    size = 1,
                    color = "orange") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank(),
             axis.title.x = element_blank()) +
       xlab(expression(scale(ln(D_IND[175]))))

#_____________________________________________________________________________________________________________
# 9c. Fall ----
#_____________________________________________________________________________________________________________

uhc.fall <- ggplot(data = sim.fall) +
       theme_bw() +
       facet_wrap(~model) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = test.fall.all[test.fall.all$case_ == 0,], 
                    aes(x = D2.indiv.hr.350.P),
                    linetype = "dashed") +
       geom_density(data = test.fall.all[test.fall.all$case_ == 1,], 
                    aes(x = D2.indiv.hr.350.P),
                    size = 1.25,
                    color = "red") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank(),
             axis.title.x = element_blank()) +
       xlab(expression(scale(ln(D_IND[350]))))

#_____________________________________________________________________________________________________________
# 9d. Dispersal ----
#_____________________________________________________________________________________________________________

uhc.disp <- ggplot(data = sim.disp) +
       theme_bw() +
       facet_wrap(~model) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = test.disp.all[test.disp.all$case_ == 0,], 
                    aes(x = D2.indiv.hr.175.P),
                    linetype = "dashed") +
       geom_density(data = test.disp.all[test.disp.all$case_ == 1,], 
                    aes(x = D2.indiv.hr.175.P),
                    size = 1,
                    color = "black") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank(),
             axis.title.x = element_blank()) +
       xlab(expression(scale(ln(D_IND[175]))))

#_____________________________________________________________________________________________________________
# 9e. All plots together ----
#_____________________________________________________________________________________________________________

cowplot::plot_grid(uhc.summer, uhc.disp,
                   nrow = 2)

#_____________________________________________________________________________________________________________
# 10. Save image ----
#_____________________________________________________________________________________________________________

save.image("uhc_plots.RData")
