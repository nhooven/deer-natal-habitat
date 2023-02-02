# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 07 - Mahalanobis distance
# Author: Nathan D. Hooven
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Email: nathan.d.hooven@gmail.com
# Date started: 5 May 2022
# Date completed: 5 May 2022
# Date modified: 6 May 2022
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data

#__________________________________________________________________________________________________
# 2. Read in data ----
#__________________________________________________________________________________________________

# pre-EHRM data
all.pre.data <- read.csv("all_pre_data_2.csv")

# create "burst_" column
all.pre.data$burst_ <- paste0(all.pre.data$Deer, ".", all.pre.data$MV.ID)

# EHRM data
ehrm.data <- read.csv("deerdata_5.csv")

#__________________________________________________________________________________________________
# 3. Calculate Mahalanobis distance of each used/available EHRM point, by EHRM ----
#__________________________________________________________________________________________________

ehrm.data.2 <- data.frame()

for (i in unique(ehrm.data$burst_)) {

    MV <- i
    
    Mv.data <- ehrm.data %>% filter(burst_ == MV)
    
    # subset frag data for each buffer size
    Mv.175 <- Mv.data %>% dplyr::select(dev.175, forest.175, open.175, 
                                        ai.175, iji.175, ed.175, prd.175)
    
    Mv.250 <- Mv.data %>% dplyr::select(dev.250, forest.250, open.250,
                                        ai.250, iji.250, ed.250, prd.250)
    
    Mv.350 <- Mv.data %>% dplyr::select(dev.350, forest.350, open.350, 
                                        ai.350, iji.350, ed.350, prd.350)
    
    Mv.500 <- Mv.data %>% dplyr::select(dev.500, forest.500, open.500, 
                                        ai.500, iji.500, ed.500, prd.500)
    
    #__________________________________________________________________________________________________
    # 3a. Distance from all deer HRs, by buffer size ----
    #__________________________________________________________________________________________________
    
    # 175-m
    all.hr.175 <- all.pre.data %>% dplyr::select(dev.175, forest.175, open.175, 
                                                 ai.175, iji.175, ed.175, prd.175)
    
    # calculate D^2 for each observation in EHRM data
    D2.all.hr.175 <- mahalanobis(x = Mv.175, center = colMeans(all.hr.175), cov = cov(all.hr.175))
    
    # 250-m
    all.hr.250 <- all.pre.data %>% dplyr::select(dev.250, forest.250, open.250, 
                                                 ai.250, iji.250, ed.250, prd.250)
    
    # calculate D^2 for each observation in EHRM data
    D2.all.hr.250 <- mahalanobis(x = Mv.250, center = colMeans(all.hr.250), cov = cov(all.hr.250))
    
    # 350-m
    all.hr.350 <- all.pre.data %>% dplyr::select(dev.350, forest.350, open.350, 
                                                 ai.350, iji.350, ed.350, prd.350)
    
    # calculate D^2 for each observation in EHRM data
    D2.all.hr.350 <- mahalanobis(x = Mv.350, center = colMeans(all.hr.350), cov = cov(all.hr.350))
    
    # 500-m
    all.hr.500 <- all.pre.data %>% dplyr::select(dev.500, forest.500, open.500, 
                                                 ai.500, iji.500, ed.500, prd.500)
    
    # calculate D^2 for each observation in EHRM data
    D2.all.hr.500 <- mahalanobis(x = Mv.500, center = colMeans(all.hr.500), cov = cov(all.hr.500))
    
    # bind together
    D2.all.hr <- data.frame(D2.all.hr.175 = D2.all.hr.175,
                            D2.all.hr.250 = D2.all.hr.250,
                            D2.all.hr.350 = D2.all.hr.350,
                            D2.all.hr.500 = D2.all.hr.500)
    
    Mv.data <- cbind(Mv.data, D2.all.hr)
    
    #__________________________________________________________________________________________________
    # 3b. Distance from individual deer HR, by buffer size ----
    #__________________________________________________________________________________________________
    
    # 175-m
    indiv.hr.175 <- all.pre.data %>% filter(burst_ == MV) %>% 
                                            dplyr::select(dev.175, forest.175, open.175, 
                                                          ai.175, iji.175, ed.175, prd.175)
    
    # calculate D^2 for each observation in EHRM data
    D2.indiv.hr.175 <- mahalanobis(x = Mv.175, center = colMeans(indiv.hr.175), cov = cov(indiv.hr.175))
    
    # 250-m
    indiv.hr.250 <- all.pre.data %>% filter(burst_ == MV) %>% 
                                     dplyr::select(dev.250, forest.250, open.250, 
                                                   ai.250, iji.250, ed.250, prd.250)
    
    # calculate D^2 for each observation in EHRM data
    D2.indiv.hr.250 <- mahalanobis(x = Mv.250, center = colMeans(indiv.hr.250), cov = cov(indiv.hr.250))
    
    # 350-m
    indiv.hr.350 <- all.pre.data %>% filter(burst_ == MV) %>% 
                                     dplyr::select(dev.350, forest.350, open.350, 
                                                   ai.350, iji.350, ed.350, prd.350)
    
    # calculate D^2 for each observation in EHRM data
    D2.indiv.hr.350 <- mahalanobis(x = Mv.350, center = colMeans(indiv.hr.350), cov = cov(indiv.hr.350))
    
    # 500-m
    indiv.hr.500 <- all.pre.data %>% filter(burst_ == MV) %>% 
                                     dplyr::select(dev.500, forest.500, open.500, 
                                                   ai.500, iji.500, ed.500, prd.500)
    
    # calculate D^2 for each observation in EHRM data
    D2.indiv.hr.500 <- mahalanobis(x = Mv.500, center = colMeans(indiv.hr.500), cov = cov(indiv.hr.500))
    
    # bind together
    D2.indiv.hr <- data.frame(D2.indiv.hr.175 = D2.indiv.hr.175,
                              D2.indiv.hr.250 = D2.indiv.hr.250,
                              D2.indiv.hr.350 = D2.indiv.hr.350,
                              D2.indiv.hr.500 = D2.indiv.hr.500)
    
    Mv.data <- cbind(Mv.data, D2.indiv.hr)
                                             
    # bind to master df
    ehrm.data.2 <- rbind(ehrm.data.2, Mv.data)
    
}

#__________________________________________________________________________________________________
# 4. Write to .csv for SSF modeling ----
#__________________________________________________________________________________________________

write.csv(ehrm.data.2, "deerdata_6.csv")
  