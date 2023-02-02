# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 09 - SSF modeling
# Author: Nathan D. Hooven
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Email: nathan.d.hooven@gmail.com
# Date started: 5 May 2022
# Date completed: 6 May 2022
# Date modified: 13 Jan 2023
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(lubridate)     # work with dates
library(survival)      # conditional logistic regressions
library(AICcmodavg)    # model selection
library(glmmTMB)       # random slope models
library(broom.mixed)   # tidy model outputs

#__________________________________________________________________________________________________
# 2. Load in data ----
#__________________________________________________________________________________________________

deer.data <- read.csv("deerdata_6.csv")

# select only columns we need for modeling, create pseudothreshold and squared variables
deer.data.1 <- deer.data %>% dplyr::select(Deer, MV.ID, burst_, step_id_, case_, sl_, ta_, inflec, t1_,
                                           dAg, dForest, dRoad, dStream,
                                           D2.all.hr.175,  D2.all.hr.250,  D2.all.hr.350,  D2.all.hr.500, 
                                           D2.indiv.hr.175,  D2.indiv.hr.250,  D2.indiv.hr.350,  D2.indiv.hr.500) %>%
                             mutate(sl_P = log(sl_ + 0.00001), 
                                    dAg.Q = dAg^2, 
                                    dForest.Q = dForest^2, 
                                    dRoad.Q = dRoad^2, 
                                    dStream.Q = dStream^2,
                                    D2.all.hr.175.Q = D2.all.hr.175^2,  
                                    D2.all.hr.250.Q = D2.all.hr.250^2,  
                                    D2.all.hr.350.Q = D2.all.hr.350^2,  
                                    D2.all.hr.500.Q = D2.all.hr.500^2, 
                                    D2.indiv.hr.175.Q = D2.indiv.hr.175^2,  
                                    D2.indiv.hr.250.Q = D2.indiv.hr.250^2,  
                                    D2.indiv.hr.350.Q = D2.indiv.hr.350^2,  
                                    D2.indiv.hr.500.Q = D2.indiv.hr.500^2,
                                    dAg.P = log(dAg + 1), 
                                    dForest.P = log(dForest + 1), 
                                    dRoad.P = log(dRoad + 1), 
                                    dStream.P = log(dStream + 1),
                                    D2.all.hr.175.P = log(D2.all.hr.175),  
                                    D2.all.hr.250.P = log(D2.all.hr.250),  
                                    D2.all.hr.350.P = log(D2.all.hr.350),  
                                    D2.all.hr.500.P = log(D2.all.hr.500), 
                                    D2.indiv.hr.175.P = log(D2.indiv.hr.175),  
                                    D2.indiv.hr.250.P = log(D2.indiv.hr.250),  
                                    D2.indiv.hr.350.P = log(D2.indiv.hr.350),  
                                    D2.indiv.hr.500.P = log(D2.indiv.hr.500))

# subset by season
deer.data.1$t1_ <- as.POSIXct(ymd_hms(deer.data.1$t1_, tz = "America/Chicago"))

# for loop to assign seasons to each EHRM
deer.data.2 <- data.frame()

for (x in unique(deer.data.1$burst_)) {
  
  burst <- x
  
  burst.data <- deer.data.1 %>% filter(burst_ == burst)
  
  # assign season
  burst.data$Season <- ifelse(substr(burst.data$t1_[1], 6, 7) %in% c("03", "04", "05") & burst.data$MV.ID != "DISP",
                              "Spring",
                              ifelse(substr(burst.data$t1_[1], 6, 7) %in% c("06", "07", "08") & burst.data$MV.ID != "DISP",
                                     "Summer",
                                     ifelse(substr(burst.data$t1_[1], 6, 7) %in% c("09", "10", "11") & burst.data$MV.ID != "DISP",
                                            "Fall",
                                            "Disp")))
  
  # bind to master df
  deer.data.2 <- rbind(deer.data.2, burst.data)
  
}

deer.data.spring <- deer.data.2 %>% filter(Season == "Spring") 
deer.data.summer <- deer.data.2 %>% filter(Season == "Summer") 
deer.data.fall <- deer.data.2 %>% filter(Season == "Fall") 
deer.data.disp <- deer.data.2 %>% filter(Season == "Disp") 

#__________________________________________________________________________________________________
# 3. Spring models ----

spring.models <- list()

#__________________________________________________________________________________________________
# 3a. Correlation and scaled variables ----
#__________________________________________________________________________________________________

cor(deer.data.spring[ ,c(10:21)])

cor(deer.data.spring[ ,c(14, 26, 42)])

write.table(cor(deer.data.spring[ ,c(10:21)]), "clipboard", sep = "\t")

deer.data.spring.s <- deer.data.spring %>% mutate_at(c(6, 10:46), scale)

#__________________________________________________________________________________________________
# 3b. dAg ----
#__________________________________________________________________________________________________

spring.models[[1]] <- clogit(case_ ~
                               dAg +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[2]] <- clogit(case_ ~
                               dAg +
                               dAg.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[3]] <- clogit(case_ ~
                               dAg.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

aictab(spring.models[1:3])

#__________________________________________________________________________________________________
# 3c. dForest ----
#__________________________________________________________________________________________________

spring.models[[4]] <- clogit(case_ ~
                               dForest +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[5]] <- clogit(case_ ~
                               dForest +
                               dForest.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[6]] <- clogit(case_ ~
                               dForest.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

aictab(spring.models[4:6])

#__________________________________________________________________________________________________
# 3d. dRoad ----
#__________________________________________________________________________________________________

spring.models[[7]] <- clogit(case_ ~
                               dRoad +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[8]] <- clogit(case_ ~
                               dRoad +
                               dRoad.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[9]] <- clogit(case_ ~
                               dRoad.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

aictab(spring.models[7:9])

#__________________________________________________________________________________________________
# 3e. dStream ----
#__________________________________________________________________________________________________

spring.models[[10]] <- clogit(case_ ~
                               dStream +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[11]] <- clogit(case_ ~
                               dStream +
                               I(dStream^2) +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[12]] <- clogit(case_ ~
                               dStream.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

aictab(spring.models[10:12])

#__________________________________________________________________________________________________
# 3f. Base model hypotheses ----
#__________________________________________________________________________________________________

spring.models[[13]] <- clogit(case_ ~
                                dAg + dAg.Q +
                                dForest.P +
                                dRoad +
                               dStream + dStream.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[14]] <- clogit(case_ ~
                                dForest.P +
                               dStream + dStream.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[15]] <- clogit(case_ ~
                               dAg + dAg.Q +
                                dRoad +
                               strata(step_id_),
                             data = deer.data.spring.s)

aictab(spring.models[13:15])

#__________________________________________________________________________________________________
# 3g. Scale/form for D2-all HRs ----
#__________________________________________________________________________________________________

spring.models[[16]] <- clogit(case_ ~
                                D2.all.hr.175 +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[17]] <- clogit(case_ ~
                                D2.all.hr.175 +
                                D2.all.hr.175.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[18]] <- clogit(case_ ~
                                D2.all.hr.175.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[19]] <- clogit(case_ ~
                                D2.all.hr.250 +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[20]] <- clogit(case_ ~
                                D2.all.hr.250 +
                                D2.all.hr.250.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[21]] <- clogit(case_ ~
                                D2.all.hr.250.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[22]] <- clogit(case_ ~
                                D2.all.hr.350 +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[23]] <- clogit(case_ ~
                                D2.all.hr.350 +
                                D2.all.hr.350.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[24]] <- clogit(case_ ~
                                D2.all.hr.350.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[25]] <- clogit(case_ ~
                                D2.all.hr.500 +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[26]] <- clogit(case_ ~
                                D2.all.hr.500 +
                                D2.all.hr.500.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[27]] <- clogit(case_ ~
                                D2.all.hr.500.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

aictab(spring.models[16:27], modnames = c("175", "175.Q", "175.P", 
                                          "250", "250.Q", "250.P",  
                                          "350", "350.Q", "350.P",
                                          "500", "500.Q", "500.P"))

write.table(aictab(spring.models[16:27], modnames = c("175", "175.Q", "175.P", 
                                          "250", "250.Q", "250.P",  
                                          "350", "350.Q", "350.P",
                                          "500", "500.Q", "500.P")),
            "clipboard",
            sep = "\t")

#__________________________________________________________________________________________________
# 3h. Scale for D2-indiv HRs ----
#__________________________________________________________________________________________________

spring.models[[28]] <- clogit(case_ ~
                                D2.indiv.hr.175 +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[29]] <- clogit(case_ ~
                                D2.indiv.hr.175 +
                                D2.indiv.hr.175.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[30]] <- clogit(case_ ~
                                D2.indiv.hr.175.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[31]] <- clogit(case_ ~
                                D2.indiv.hr.250 +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[32]] <- clogit(case_ ~
                                D2.indiv.hr.250 +
                                D2.indiv.hr.250.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[33]] <- clogit(case_ ~
                                D2.indiv.hr.250.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[34]] <- clogit(case_ ~
                                D2.indiv.hr.350 +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[35]] <- clogit(case_ ~
                                D2.indiv.hr.350 +
                                D2.indiv.hr.350.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[36]] <- clogit(case_ ~
                                D2.indiv.hr.350.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[37]] <- clogit(case_ ~
                                D2.indiv.hr.500 +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[38]] <- clogit(case_ ~
                                D2.indiv.hr.500 +
                                D2.indiv.hr.500.Q +
                               strata(step_id_),
                             data = deer.data.spring.s)

spring.models[[39]] <- clogit(case_ ~
                                D2.indiv.hr.500.P +
                               strata(step_id_),
                             data = deer.data.spring.s)

write.table(aictab(spring.models[28:39], modnames = c("175", "175.Q", "175.P", 
                                          "250", "250.Q", "250.P",  
                                          "350", "350.Q", "350.P",
                                          "500", "500.Q", "500.P")),
            "clipboard",
            sep = "\t")

#__________________________________________________________________________________________________
# 3i. Final model selection ----
#__________________________________________________________________________________________________

spring.models.40.struc <- glmmTMB(case_ ~
                                dForest.P +
                                dStream + dStream.Q +
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  dForest.P | Deer) +
                                  (0 +  dStream | Deer),
                           family = poisson,
                             data = deer.data.spring.s,
                           doFit = FALSE)

spring.models.40.struc$parameters$theta[1] = log(1e3)
spring.models.40.struc$mapArg = list(theta = factor(c(NA, 1:2)))

spring.models[[40]] <- glmmTMB:::fitTMB(spring.models.40.struc)

spring.models.41.struc <- glmmTMB(case_ ~
                                dForest.P +
                                dStream + dStream.Q +
                                D2.all.hr.175 +
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  dForest.P | Deer) +
                                  (0 +  dStream | Deer) +
                                  (0 +  D2.all.hr.175 | Deer),
                           family = poisson,
                             data = deer.data.spring.s,
                           doFit = FALSE)

spring.models.41.struc$parameters$theta[1] = log(1e3)
spring.models.41.struc$mapArg = list(theta = factor(c(NA, 1:3)))

spring.models[[41]] <- glmmTMB:::fitTMB(spring.models.41.struc)

spring.models.42.struc <- glmmTMB(case_ ~
                                dForest.P +
                                dStream + dStream.Q +
                                D2.all.hr.175 +
                                D2.indiv.hr.175.P + 
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  dForest.P | Deer) +
                                  (0 +  dStream | Deer) +
                                  (0 +  D2.all.hr.175 | Deer) +
                                (0 +  D2.indiv.hr.175.P | Deer),
                           family = poisson,
                             data = deer.data.spring.s,
                           doFit = FALSE)

spring.models.42.struc$parameters$theta[1] = log(1e3)
spring.models.42.struc$mapArg = list(theta = factor(c(NA, 1:4)))

spring.models[[42]] <- glmmTMB:::fitTMB(spring.models.42.struc)

aictab(spring.models[40:42], modnames = c("base", "base + all", "base + all + indiv"))

summary(spring.models[[42]])

write.table(cor(deer.data.spring.s[, c("dForest.P",
                           "dStream",
                           "dStream.Q",
                           "D2.all.hr.175", 
                           "D2.indiv.hr.175.P")]),
            "clipboard",
            sep = "\t")

# manual VIF
# dForest.P
1 / (1 - summary(lm(dForest.P ~ dStream + D2.all.hr.175 + D2.indiv.hr.175.P, data = deer.data.spring.s))$r.squared)

# dStream
1 / (1 - summary(lm(dStream ~ dForest.P + D2.all.hr.175 + D2.indiv.hr.175.P, data = deer.data.spring.s))$r.squared)

# D2.all.hr.175.P
1 / (1 - summary(lm(D2.all.hr.175.P ~ dForest.P + dStream + D2.indiv.hr.175.P, data = deer.data.spring.s))$r.squared)

# D2.indiv.hr.175.P
1 / (1 - summary(lm(D2.indiv.hr.175.P ~ dForest.P + dStream + D2.all.hr.175.P, data = deer.data.spring.s))$r.squared)


ggplot(data = deer.data.spring.s, aes(D2.indiv.hr.500, color = case_)) +
  geom_density() +
  coord_cartesian(xlim = c(-0.5, 2))

#__________________________________________________________________________________________________
# 4. Summer models ----

summer.models <- list()

#__________________________________________________________________________________________________
# 4a. Correlation and scaled variables ----
#__________________________________________________________________________________________________

write.table(cor(deer.data.summer[ ,c(10:21)]),
            "clipboard",
            sep = "\t")

deer.data.summer.s <- deer.data.summer %>% mutate_at(c(6, 10:46), scale)

#__________________________________________________________________________________________________
# 4b. dAg ----
#__________________________________________________________________________________________________

summer.models[[1]] <- clogit(case_ ~
                               dAg +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[2]] <- clogit(case_ ~
                               dAg +
                               dAg.Q +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[3]] <- clogit(case_ ~
                               dAg.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

aictab(summer.models[1:3])

#__________________________________________________________________________________________________
# 4c. dForest ----
#__________________________________________________________________________________________________

summer.models[[4]] <- clogit(case_ ~
                               dForest +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[5]] <- clogit(case_ ~
                               dForest +
                               dForest.Q +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[6]] <- clogit(case_ ~
                               dForest.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

aictab(summer.models[4:6])

#__________________________________________________________________________________________________
# 4d. dRoad ----
#__________________________________________________________________________________________________

summer.models[[7]] <- clogit(case_ ~
                               dRoad +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[8]] <- clogit(case_ ~
                               dRoad +
                               dRoad.Q +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[9]] <- clogit(case_ ~
                               dRoad.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

aictab(summer.models[7:9])

#__________________________________________________________________________________________________
# 4e. dStream ----
#__________________________________________________________________________________________________

summer.models[[10]] <- clogit(case_ ~
                               dStream +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[11]] <- clogit(case_ ~
                               dStream +
                               I(dStream^2) +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[12]] <- clogit(case_ ~
                               dStream.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

aictab(summer.models[10:12])

#__________________________________________________________________________________________________
# 4f. Base model hypotheses ----
#__________________________________________________________________________________________________

summer.models[[13]] <- clogit(case_ ~
                                dAg.P +
                                dForest.P +
                                dRoad.P +
                               dStream +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[14]] <- clogit(case_ ~
                                dForest.P +
                               dStream +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[15]] <- clogit(case_ ~
                               dAg.P +
                                dRoad.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

aictab(summer.models[13:15])

#__________________________________________________________________________________________________
# 4g. Scale/form for D2-all HRs ----
#__________________________________________________________________________________________________

summer.models[[16]] <- clogit(case_ ~
                                D2.all.hr.175 +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[17]] <- clogit(case_ ~
                                D2.all.hr.175 +
                                D2.all.hr.175.Q +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[18]] <- clogit(case_ ~
                                D2.all.hr.175.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[19]] <- clogit(case_ ~
                                D2.all.hr.250 +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[20]] <- clogit(case_ ~
                                D2.all.hr.250 +
                                D2.all.hr.250.Q +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[21]] <- clogit(case_ ~
                                D2.all.hr.250.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[22]] <- clogit(case_ ~
                                D2.all.hr.350 +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[23]] <- clogit(case_ ~
                                D2.all.hr.350 +
                                D2.all.hr.350.Q +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[24]] <- clogit(case_ ~
                                D2.all.hr.350.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[25]] <- clogit(case_ ~
                                D2.all.hr.500 +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[26]] <- clogit(case_ ~
                                D2.all.hr.500 +
                                D2.all.hr.500.Q +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[27]] <- clogit(case_ ~
                                D2.all.hr.500.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

write.table(aictab(summer.models[16:27], modnames = c("175", "175.Q", "175.P", 
                                          "250", "250.Q", "250.P",  
                                          "350", "350.Q", "350.P",
                                          "500", "500.Q", "500.P")),
            "clipboard",
            sep = "\t")

#__________________________________________________________________________________________________
# 4h. Scale for D2-indiv HRs ----
#__________________________________________________________________________________________________

summer.models[[28]] <- clogit(case_ ~
                                D2.indiv.hr.175 +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[29]] <- clogit(case_ ~
                                D2.indiv.hr.175 +
                                D2.indiv.hr.175.Q +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[30]] <- clogit(case_ ~
                                D2.indiv.hr.175.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[31]] <- clogit(case_ ~
                                D2.indiv.hr.250 +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[32]] <- clogit(case_ ~
                                D2.indiv.hr.250 +
                                D2.indiv.hr.250.Q +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[33]] <- clogit(case_ ~
                                D2.indiv.hr.250.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[34]] <- clogit(case_ ~
                                D2.indiv.hr.350 +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[35]] <- clogit(case_ ~
                                D2.indiv.hr.350 +
                                D2.indiv.hr.350.Q +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[36]] <- clogit(case_ ~
                                D2.indiv.hr.350.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[37]] <- clogit(case_ ~
                                D2.indiv.hr.500 +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[38]] <- clogit(case_ ~
                                D2.indiv.hr.500 +
                                D2.indiv.hr.500.Q +
                               strata(step_id_),
                             data = deer.data.summer.s)

summer.models[[39]] <- clogit(case_ ~
                                D2.indiv.hr.500.P +
                               strata(step_id_),
                             data = deer.data.summer.s)

write.table(aictab(summer.models[28:39], modnames = c("175", "175.Q", "175.P", 
                                          "250", "250.Q", "250.P",  
                                          "350", "350.Q", "350.P",
                                          "500", "500.Q", "500.P")),
            "clipboard",
            sep = "\t")

#__________________________________________________________________________________________________
# 4i. Final model selection ----
#__________________________________________________________________________________________________

summer.models.40.struc <- glmmTMB(case_ ~
                                dForest.P +
                                dStream + 
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  dForest.P | Deer) +
                                  (0 +  dStream | Deer),
                           family = poisson,
                             data = deer.data.summer.s,
                           doFit = FALSE)

summer.models.40.struc$parameters$theta[1] = log(1e3)
summer.models.40.struc$mapArg = list(theta = factor(c(NA, 1:2)))

summer.models[[40]] <- glmmTMB:::fitTMB(summer.models.40.struc)

summer.models.41.struc <- glmmTMB(case_ ~
                                dForest.P +
                                dStream + 
                                D2.all.hr.250.P +
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  dForest.P | Deer) +
                                  (0 +  dStream | Deer) +
                                  (0 +  D2.all.hr.250.P | Deer),
                           family = poisson,
                             data = deer.data.summer.s,
                           doFit = FALSE)

summer.models.41.struc$parameters$theta[1] = log(1e3)
summer.models.41.struc$mapArg = list(theta = factor(c(NA, 1:3)))

summer.models[[41]] <- glmmTMB:::fitTMB(summer.models.41.struc)

summer.models.42.struc <- glmmTMB(case_ ~
                                dForest.P +
                                dStream + 
                                D2.all.hr.250.P +
                                D2.indiv.hr.175.P + 
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  dForest.P | Deer) +
                                  (0 +  dStream | Deer) +
                                  (0 +  D2.all.hr.250.P | Deer) +
                                  (0 +  D2.indiv.hr.175.P | Deer),
                           family = poisson,
                             data = deer.data.summer.s,
                           doFit = FALSE)

summer.models.42.struc$parameters$theta[1] = log(1e3)
summer.models.42.struc$mapArg = list(theta = factor(c(NA, 1:4)))

summer.models[[42]] <- glmmTMB:::fitTMB(summer.models.42.struc)

aictab(summer.models[40:42], modnames = c("base", "base + all", "base + all + indiv"))

summary(summer.models[[42]])

write.table(cor(deer.data.spring.s[, c("dForest.P",
                                       "dStream",
                                       "D2.all.hr.250.P", 
                                       "D2.indiv.hr.175.P")]),
            "clipboard",
            sep = "\t")

# manual VIF
# dForest.P
1 / (1 - summary(lm(dForest.P ~ dStream + D2.all.hr.250.P + D2.indiv.hr.175.P, data = deer.data.summer.s))$r.squared)

# dStream
1 / (1 - summary(lm(dStream ~ dForest.P + D2.all.hr.250.P + D2.indiv.hr.175.P, data = deer.data.spring.s))$r.squared)

# D2.all.hr.250.P
1 / (1 - summary(lm(D2.all.hr.250.P ~ dForest.P + dStream + D2.indiv.hr.175.P, data = deer.data.spring.s))$r.squared)

# D2.indiv.hr.175.P
1 / (1 - summary(lm(D2.indiv.hr.175.P ~ dForest.P + dStream + D2.all.hr.250.P, data = deer.data.spring.s))$r.squared)

ggplot(data = deer.data.summer.s, aes(D2.indiv.hr.175.P, color = case_)) +
  geom_density()

#__________________________________________________________________________________________________
# 5. Fall models ----

fall.models <- list()

#__________________________________________________________________________________________________
# 5a. Correlation and scaled variables ----
#__________________________________________________________________________________________________

write.table(cor(deer.data.fall[ ,c(10:21)]),
            "clipboard",
            sep = "\t")

deer.data.fall.s <- deer.data.fall %>% mutate_at(c(6, 10:46), scale)

#__________________________________________________________________________________________________
# 5b. dAg ----
#__________________________________________________________________________________________________

fall.models[[1]] <- clogit(case_ ~
                               dAg +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[2]] <- clogit(case_ ~
                               dAg +
                               dAg.Q +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[3]] <- clogit(case_ ~
                               dAg.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

aictab(fall.models[1:3])

#__________________________________________________________________________________________________
# 5c. dForest ----
#__________________________________________________________________________________________________

fall.models[[4]] <- clogit(case_ ~
                               dForest +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[5]] <- clogit(case_ ~
                               dForest +
                               dForest.Q +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[6]] <- clogit(case_ ~
                               dForest.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

aictab(fall.models[4:6])

#__________________________________________________________________________________________________
# 5d. dRoad ----
#__________________________________________________________________________________________________

fall.models[[7]] <- clogit(case_ ~
                               dRoad +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[8]] <- clogit(case_ ~
                               dRoad +
                               dRoad.Q +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[9]] <- clogit(case_ ~
                               dRoad.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

aictab(fall.models[7:9])

#__________________________________________________________________________________________________
# 5e. dStream ----
#__________________________________________________________________________________________________

fall.models[[10]] <- clogit(case_ ~
                               dStream +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[11]] <- clogit(case_ ~
                               dStream +
                               I(dStream^2) +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[12]] <- clogit(case_ ~
                               dStream.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

aictab(fall.models[10:12])

#__________________________________________________________________________________________________
# 5f. Base model hypotheses ----
#__________________________________________________________________________________________________

fall.models[[13]] <- clogit(case_ ~
                                dAg.P +
                                dForest.P +
                                dRoad.P +
                               dStream +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[14]] <- clogit(case_ ~
                                dForest.P +
                               dStream +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[15]] <- clogit(case_ ~
                               dAg.P +
                                dRoad.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

aictab(fall.models[13:15])

#__________________________________________________________________________________________________
# 5g. Scale/form for D2-all HRs ----
#__________________________________________________________________________________________________

fall.models[[16]] <- clogit(case_ ~
                                D2.all.hr.175 +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[17]] <- clogit(case_ ~
                                D2.all.hr.175 +
                                D2.all.hr.175.Q +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[18]] <- clogit(case_ ~
                                D2.all.hr.175.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[19]] <- clogit(case_ ~
                                D2.all.hr.250 +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[20]] <- clogit(case_ ~
                                D2.all.hr.250 +
                                D2.all.hr.250.Q +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[21]] <- clogit(case_ ~
                                D2.all.hr.250.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[22]] <- clogit(case_ ~
                                D2.all.hr.350 +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[23]] <- clogit(case_ ~
                                D2.all.hr.350 +
                                D2.all.hr.350.Q +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[24]] <- clogit(case_ ~
                                D2.all.hr.350.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[25]] <- clogit(case_ ~
                                D2.all.hr.500 +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[26]] <- clogit(case_ ~
                                D2.all.hr.500 +
                                D2.all.hr.500.Q +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[27]] <- clogit(case_ ~
                                D2.all.hr.500.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

write.table(aictab(fall.models[16:27], modnames = c("175", "175.Q", "175.P", 
                                          "250", "250.Q", "250.P",  
                                          "350", "350.Q", "350.P",
                                          "500", "500.Q", "500.P")),
            "clipboard",
            sep = "\t")

#__________________________________________________________________________________________________
# 5h. Scale for D2-indiv HRs ----
#__________________________________________________________________________________________________

fall.models[[28]] <- clogit(case_ ~
                                D2.indiv.hr.175 +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[29]] <- clogit(case_ ~
                                D2.indiv.hr.175 +
                                D2.indiv.hr.175.Q +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[30]] <- clogit(case_ ~
                                D2.indiv.hr.175.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[31]] <- clogit(case_ ~
                                D2.indiv.hr.250 +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[32]] <- clogit(case_ ~
                                D2.indiv.hr.250 +
                                D2.indiv.hr.250.Q +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[33]] <- clogit(case_ ~
                                D2.indiv.hr.250.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[34]] <- clogit(case_ ~
                                D2.indiv.hr.350 +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[35]] <- clogit(case_ ~
                                D2.indiv.hr.350 +
                                D2.indiv.hr.350.Q +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[36]] <- clogit(case_ ~
                                D2.indiv.hr.350.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[37]] <- clogit(case_ ~
                                D2.indiv.hr.500 +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[38]] <- clogit(case_ ~
                                D2.indiv.hr.500 +
                                D2.indiv.hr.500.Q +
                               strata(step_id_),
                             data = deer.data.fall.s)

fall.models[[39]] <- clogit(case_ ~
                                D2.indiv.hr.500.P +
                               strata(step_id_),
                             data = deer.data.fall.s)

write.table(aictab(fall.models[28:39], modnames = c("175", "175.Q", "175.P", 
                                          "250", "250.Q", "250.P",  
                                          "350", "350.Q", "350.P",
                                          "500", "500.Q", "500.P")),
            "clipboard",
            sep = "\t")

#__________________________________________________________________________________________________
# 5i. Final model selection ----
#__________________________________________________________________________________________________

fall.models.40.struc <- glmmTMB(case_ ~
                                dAg.P +
                                dForest.P +
                                dRoad.P +
                                dStream +
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  dForest.P | Deer) +
                                  (0 +  dStream | Deer),
                           family = poisson,
                             data = deer.data.fall.s,
                           doFit = FALSE)

fall.models.40.struc$parameters$theta[1] = log(1e3)
fall.models.40.struc$mapArg = list(theta = factor(c(NA, 1:2)))

fall.models[[40]] <- glmmTMB:::fitTMB(fall.models.40.struc)

fall.models.41.struc <- glmmTMB(case_ ~
                                dAg.P +
                                dForest.P +
                                dRoad.P +
                                dStream +
                                  D2.all.hr.175.P +
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  dForest.P | Deer) +
                                  (0 +  dStream | Deer) +
                                  (0 +  D2.all.hr.175.P | Deer),
                           family = poisson,
                             data = deer.data.fall.s,
                           doFit = FALSE)

fall.models.41.struc$parameters$theta[1] = log(1e3)
fall.models.41.struc$mapArg = list(theta = factor(c(NA, 1:3)))

fall.models[[41]] <- glmmTMB:::fitTMB(fall.models.41.struc)

fall.models.42.struc <- glmmTMB(case_ ~
                                dAg.P +
                                dForest.P +
                                dRoad.P +
                                dStream +
                                  D2.all.hr.175.P +
                                  D2.indiv.hr.350.P +
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  dForest.P | Deer) +
                                  (0 +  dStream | Deer) +
                                  (0 +  D2.indiv.hr.350.P | Deer),
                           family = poisson,
                             data = deer.data.fall.s,
                           doFit = FALSE)

fall.models.42.struc$parameters$theta[1] = log(1e3)
fall.models.42.struc$mapArg = list(theta = factor(c(NA, 1:3)))

fall.models[[42]] <- glmmTMB:::fitTMB(fall.models.42.struc)

aictab(fall.models[40:42], modnames = c("base", "base + all", "base + all + indiv"))

summary(fall.models[[42]])

write.table(cor(deer.data.spring.s[, c("dAg.P",
                                       "dForest.P",
                                       "dRoad.P",
                                       "dStream",
                                       "D2.all.hr.175.P", 
                                       "D2.indiv.hr.350.P")]),
            "clipboard",
            sep = "\t")

# manual VIF
# dForest.P
1 / (1 - summary(lm(dAg.P ~ dForest.P + dRoad.P + dStream + D2.all.hr.175.P + D2.indiv.hr.350.P, data = deer.data.fall.s))$r.squared)

# dForest.P
1 / (1 - summary(lm(dForest.P ~ dAg.P + dRoad.P + dStream + D2.all.hr.175.P + D2.indiv.hr.350.P, data = deer.data.fall.s))$r.squared)

# dStream
1 / (1 - summary(lm(dStream ~ dAg.P + dForest.P + dRoad.P  + D2.all.hr.175.P + D2.indiv.hr.350.P, data = deer.data.fall.s))$r.squared)

# D2.all.hr.175.P
1 / (1 - summary(lm(D2.all.hr.175.P ~ dAg.P + dForest.P + dRoad.P + dStream + D2.indiv.hr.350.P, data = deer.data.fall.s))$r.squared)

# D2.indiv.hr.350.P
1 / (1 - summary(lm(D2.indiv.hr.350.P ~ dAg.P + dForest.P + dRoad.P + dStream + D2.all.hr.175.P, data = deer.data.fall.s))$r.squared)

ggplot(data = deer.data.fall.s, aes(D2.indiv.hr.350.P, color = case_)) +
  geom_density()

#__________________________________________________________________________________________________
# 6. Dispersal models ----

disp.models <- list()

#__________________________________________________________________________________________________
# 6a. Correlation and scaled variables ----
#__________________________________________________________________________________________________

write.table(cor(deer.data.disp[ ,c(10:21)]),
            "clipboard",
            sep = "\t")

deer.data.disp.s <- deer.data.disp %>% mutate_at(c(6, 10:46), scale)

#__________________________________________________________________________________________________
# 6b. dAg ----
#__________________________________________________________________________________________________

disp.models[[1]] <- clogit(case_ ~
                               dAg +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[2]] <- clogit(case_ ~
                               dAg +
                               dAg.Q +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[3]] <- clogit(case_ ~
                               dAg.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

aictab(disp.models[1:3])

#__________________________________________________________________________________________________
# 6c. dForest ----
#__________________________________________________________________________________________________

disp.models[[4]] <- clogit(case_ ~
                               dForest +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[5]] <- clogit(case_ ~
                               dForest +
                               dForest.Q +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[6]] <- clogit(case_ ~
                               dForest.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

aictab(disp.models[4:6])

#__________________________________________________________________________________________________
# 6d. dRoad ----
#__________________________________________________________________________________________________

disp.models[[7]] <- clogit(case_ ~
                               dRoad +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[8]] <- clogit(case_ ~
                               dRoad +
                               dRoad.Q +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[9]] <- clogit(case_ ~
                               dRoad.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

aictab(disp.models[7:9])

#__________________________________________________________________________________________________
# 6e. dStream ----
#__________________________________________________________________________________________________

disp.models[[10]] <- clogit(case_ ~
                               dStream +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[11]] <- clogit(case_ ~
                               dStream +
                               I(dStream^2) +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[12]] <- clogit(case_ ~
                               dStream.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

aictab(disp.models[10:12])

#__________________________________________________________________________________________________
# 6f. Base model hypotheses ----
#__________________________________________________________________________________________________

disp.models[[13]] <- clogit(case_ ~
                                dAg.P  +
                                dForest.P +
                                dRoad.P +
                               dStream +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[14]] <- clogit(case_ ~
                                dForest.P +
                               dStream +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[15]] <- clogit(case_ ~
                               dAg.P +
                                dRoad.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

aictab(disp.models[13:15])

#__________________________________________________________________________________________________
# 6g. Scale/form for D2-all HRs ----
#__________________________________________________________________________________________________

disp.models[[16]] <- clogit(case_ ~
                                D2.all.hr.175 +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[17]] <- clogit(case_ ~
                                D2.all.hr.175 +
                                D2.all.hr.175.Q +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[18]] <- clogit(case_ ~
                                D2.all.hr.175.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[19]] <- clogit(case_ ~
                                D2.all.hr.250 +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[20]] <- clogit(case_ ~
                                D2.all.hr.250 +
                                D2.all.hr.250.Q +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[21]] <- clogit(case_ ~
                                D2.all.hr.250.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[22]] <- clogit(case_ ~
                                D2.all.hr.350 +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[23]] <- clogit(case_ ~
                                D2.all.hr.350 +
                                D2.all.hr.350.Q +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[24]] <- clogit(case_ ~
                                D2.all.hr.350.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[25]] <- clogit(case_ ~
                                D2.all.hr.500 +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[26]] <- clogit(case_ ~
                                D2.all.hr.500 +
                                D2.all.hr.500.Q +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[27]] <- clogit(case_ ~
                                D2.all.hr.500.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

write.table(aictab(disp.models[16:27], modnames = c("175", "175.Q", "175.P", 
                                          "250", "250.Q", "250.P",  
                                          "350", "350.Q", "350.P",
                                          "500", "500.Q", "500.P")),
            "clipboard",
            sep = "\t")

#__________________________________________________________________________________________________
# 6h. Scale for D2-indiv HRs ----
#__________________________________________________________________________________________________

disp.models[[28]] <- clogit(case_ ~
                                D2.indiv.hr.175 +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[29]] <- clogit(case_ ~
                                D2.indiv.hr.175 +
                                D2.indiv.hr.175.Q +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[30]] <- clogit(case_ ~
                                D2.indiv.hr.175.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[31]] <- clogit(case_ ~
                                D2.indiv.hr.250 +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[32]] <- clogit(case_ ~
                                D2.indiv.hr.250 +
                                D2.indiv.hr.250.Q +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[33]] <- clogit(case_ ~
                                D2.indiv.hr.250.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[34]] <- clogit(case_ ~
                                D2.indiv.hr.350 +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[35]] <- clogit(case_ ~
                                D2.indiv.hr.350 +
                                D2.indiv.hr.350.Q +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[36]] <- clogit(case_ ~
                                D2.indiv.hr.350.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[37]] <- clogit(case_ ~
                                D2.indiv.hr.500 +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[38]] <- clogit(case_ ~
                                D2.indiv.hr.500 +
                                D2.indiv.hr.500.Q +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[39]] <- clogit(case_ ~
                                D2.indiv.hr.500.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

write.table(aictab(disp.models[28:39], modnames = c("175", "175.Q", "175.P", 
                                          "250", "250.Q", "250.P",  
                                          "350", "350.Q", "350.P",
                                          "500", "500.Q", "500.P")),
            "clipboard",
            sep = "\t")

summary(disp.models[[30]])

#__________________________________________________________________________________________________
# 6i. Final model selection ----
#__________________________________________________________________________________________________

disp.models[[40]] <- clogit(case_ ~
                               dForest.P +
                               dStream.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[41]] <- clogit(case_ ~
                               dForest.P +
                               dStream.P +
                               D2.all.hr.175.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models[[42]] <- clogit(case_ ~
                               dForest.P +
                               dStream.P +
                               D2.all.hr.175.P +
                               D2.indiv.hr.175.P +
                               strata(step_id_),
                             data = deer.data.disp.s)

disp.models.40.struc <- glmmTMB(case_ ~
                                dForest.P +
                                dStream.P +
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  dForest.P | Deer) +
                                  (0 +  dStream.P | Deer),
                           family = poisson,
                             data = deer.data.disp.s,
                           doFit = FALSE)

disp.models.40.struc$parameters$theta[1] = log(1e3)
disp.models.40.struc$mapArg = list(theta = factor(c(NA, 1:2)))

disp.models[[40]] <- glmmTMB:::fitTMB(disp.models.40.struc)

disp.models.41.struc <- glmmTMB(case_ ~
                                dForest.P +
                                dStream.P +
                                D2.all.hr.175 + 
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  dForest.P | Deer) +
                                  (0 +  dStream.P | Deer) +
                                  (0 +  D2.all.hr.175 | Deer),
                           family = poisson,
                             data = deer.data.disp.s,
                           doFit = FALSE)

disp.models.41.struc$parameters$theta[1] = log(1e3)
disp.models.41.struc$mapArg = list(theta = factor(c(NA, 1:3)))

disp.models[[41]] <- glmmTMB:::fitTMB(disp.models.41.struc)

disp.models.42.struc <- glmmTMB(case_ ~
                                dForest.P +
                                dStream.P +
                                D2.all.hr.175 +
                                  D2.indiv.hr.175.P +
                                 sl_ + 
                                (1 | step_id_) +
                                  (0 +  D2.indiv.hr.175.P | Deer),
                           family = poisson,
                             data = deer.data.disp.s,
                           doFit = FALSE)

disp.models.42.struc$parameters$theta[1] = log(1e3)
disp.models.42.struc$mapArg = list(theta = factor(c(NA, 1)))

disp.models[[42]] <- glmmTMB:::fitTMB(disp.models.42.struc)

aictab(disp.models[40:42], modnames = c("base", "base + all", "base + all + indiv"))

summary(disp.models[[42]])

performance::check_collinearity(disp.models[[42]])

write.table(cor(deer.data.spring.s[, c("dForest.P",
                                       "dStream.P",
                                       "D2.all.hr.175", 
                                       "D2.indiv.hr.175.P")]),
            "clipboard",
            sep = "\t")

# manual VIF
# dForest.P
1 / (1 - summary(lm(dForest.P ~ dStream.P + D2.all.hr.175 + D2.indiv.hr.175.P, data = deer.data.disp.s))$r.squared)

# dStream.P
1 / (1 - summary(lm(dStream.P ~ dForest.P + D2.all.hr.175 + D2.indiv.hr.175.P, data = deer.data.disp.s))$r.squared)

# D2.all.hr.175.P
1 / (1 - summary(lm(D2.all.hr.175 ~ dForest.P + dStream + D2.indiv.hr.175.P, data = deer.data.disp.s))$r.squared)

# D2.indiv.hr.175.P
1 / (1 - summary(lm(D2.indiv.hr.175.P ~ dForest.P + dStream + D2.all.hr.175, data = deer.data.disp.s))$r.squared)

ggplot(data = deer.data.disp.s, aes(D2.indiv.hr.175.P, color = case_)) +
  geom_density()

#__________________________________________________________________________________________________
# 7. Write AIC tables ----
#__________________________________________________________________________________________________

write.table(aictab(spring.models[40:42], modnames = c("base", "base + all", "base + all + indiv")),
            "clipboard", sep = "\t")

write.table(aictab(summer.models[40:42], modnames = c("base", "base + all", "base + all + indiv")),
            "clipboard", sep = "\t")

write.table(aictab(fall.models[40:42], modnames = c("base", "base + all", "base + all + indiv")),
            "clipboard", sep = "\t")

write.table(aictab(disp.models[40:42], modnames = c("base", "base + all", "base + all + indiv")),
            "clipboard", sep = "\t")

#__________________________________________________________________________________________________
# 8. Table of parameter estimates and CIs ----
#__________________________________________________________________________________________________

# spring
spring.params <- data.frame()

for (x in 40:42) {
  
  spring.tidy <- tidy(spring.models[[x]])
  
  spring.tidy <- spring.tidy %>% mutate(season = "Spring") %>%
                                 dplyr::filter(effect == "fixed",
                                               term != "(Intercept)") %>%
                                 select(season, term, estimate, std.error) %>%
                                 mutate(CI.low = estimate - std.error*1.96,
                                        CI.upp = estimate + std.error*1.96) %>%
                                 mutate(est.CI = paste0(round(estimate, digits = 3),
                                                        " (",
                                                        round(CI.low, digits = 3), 
                                                        ", ", 
                                                        round(CI.upp, digits = 3),
                                                        ")"))
  
  spring.tidy.pivot <- spring.tidy %>% dplyr::select(term, est.CI) %>%
                                              pivot_wider(names_from = term,
                                                          values_from = est.CI) %>%
                                              mutate(model = ifelse(x == 40,
                                                                    "base",
                                                                    ifelse(x == 41,
                                                                           "base + all",
                                                                           "base + all + indiv")))
  
  spring.params <- bind_rows(spring.params, spring.tidy.pivot)
  
}

# summer
summer.params <- data.frame()

for (x in 40:42) {
  
  summer.tidy <- tidy(summer.models[[x]])
  
  summer.tidy <- summer.tidy %>% mutate(season = "summer") %>%
                                 dplyr::filter(effect == "fixed",
                                               term != "(Intercept)") %>%
                                 select(season, term, estimate, std.error) %>%
                                 mutate(CI.low = estimate - std.error*1.96,
                                        CI.upp = estimate + std.error*1.96) %>%
                                 mutate(est.CI = paste0(round(estimate, digits = 3),
                                                        " (",
                                                        round(CI.low, digits = 3), 
                                                        ", ", 
                                                        round(CI.upp, digits = 3),
                                                        ")"))
  
  summer.tidy.pivot <- summer.tidy %>% dplyr::select(term, est.CI) %>%
                                              pivot_wider(names_from = term,
                                                          values_from = est.CI) %>%
                                              mutate(model = ifelse(x == 40,
                                                                    "base",
                                                                    ifelse(x == 41,
                                                                           "base + all",
                                                                           "base + all + indiv")))
  
  summer.params <- bind_rows(summer.params, summer.tidy.pivot)
  
}

# fall
fall.params <- data.frame()

for (x in 40:42) {
  
  fall.tidy <- tidy(fall.models[[x]])
  
  fall.tidy <- fall.tidy %>% mutate(season = "fall") %>%
                                 dplyr::filter(effect == "fixed",
                                               term != "(Intercept)") %>%
                                 select(season, term, estimate, std.error) %>%
                                 mutate(CI.low = estimate - std.error*1.96,
                                        CI.upp = estimate + std.error*1.96) %>%
                                 mutate(est.CI = paste0(round(estimate, digits = 3),
                                                        " (",
                                                        round(CI.low, digits = 3), 
                                                        ", ", 
                                                        round(CI.upp, digits = 3),
                                                        ")"))
  
  fall.tidy.pivot <- fall.tidy %>% dplyr::select(term, est.CI) %>%
                                              pivot_wider(names_from = term,
                                                          values_from = est.CI) %>%
                                              mutate(model = ifelse(x == 40,
                                                                    "base",
                                                                    ifelse(x == 41,
                                                                           "base + all",
                                                                           "base + all + indiv")))
  
  fall.params <- bind_rows(fall.params, fall.tidy.pivot)
  
}

# disp
disp.params <- data.frame()

for (x in 40:42) {
  
  disp.tidy <- tidy(disp.models[[x]])
  
  disp.tidy <- disp.tidy %>% mutate(season = "disp") %>%
                                 dplyr::filter(effect == "fixed",
                                               term != "(Intercept)") %>%
                                 select(season, term, estimate, std.error) %>%
                                 mutate(CI.low = estimate - std.error*1.96,
                                        CI.upp = estimate + std.error*1.96) %>%
                                 mutate(est.CI = paste0(round(estimate, digits = 3),
                                                        " (",
                                                        round(CI.low, digits = 3), 
                                                        ", ", 
                                                        round(CI.upp, digits = 3),
                                                        ")"))
  
  disp.tidy.pivot <- disp.tidy %>% dplyr::select(term, est.CI) %>%
                                              pivot_wider(names_from = term,
                                                          values_from = est.CI) %>%
                                              mutate(model = ifelse(x == 40,
                                                                    "base",
                                                                    ifelse(x == 41,
                                                                           "base + all",
                                                                           "base + all + indiv")))
  
  disp.params <- bind_rows(disp.params, disp.tidy.pivot)
  
}

# bind all together
all.params <- bind_rows(spring.params, summer.params, fall.params, disp.params)

write.table(all.params, "clipboard", sep = "\t")

#__________________________________________________________________________________________________
# 8. Save image ----
#__________________________________________________________________________________________________

save.image(file = "SSF_models_NEW.RData")

