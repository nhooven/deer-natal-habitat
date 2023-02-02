# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 05 - Identify inflection points
# Author: Nathan D. Hooven
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Email: nathan.d.hooven@gmail.com
# Date started: 9 Dec 2020
# Date completed: 9 Dec 2020
# Date modified: 5 May 2022
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(sp)            # work with spatial data
library(raster)        # distances

#__________________________________________________________________________________________________
# 2. Read in data ----
#__________________________________________________________________________________________________

all.data <- read.csv("deerdata_3.csv")

#__________________________________________________________________________________________________
# 3. Read in data ----
#__________________________________________________________________________________________________

all.data.inflec <- data.frame()

for (i in unique(all.data$Deer)) {
  
  DeerID <- i
  
  # subset data
  Deer.data <- all.data %>% filter(Deer == DeerID)
  
  for (j in unique(Deer.data$MV.ID)) {
    
    MV <- j
    
    Mv.data <- Deer.data %>% filter(MV.ID == MV)
    
    dists <- pointDistance(Mv.data[, c("x", "y")], lonlat = FALSE)
    
    Mv.data$dists <- dists[,1]
    
    print(ggplot(data = Mv.data[Mv.data$case_ == TRUE,], aes(x = step_id_, y = dists)) +
          geom_point() +
          geom_hline(yintercept = max(Mv.data$dists[Mv.data$case_ == TRUE])) +
          theme_bw() +
          ggtitle(paste0(DeerID, ".", MV)))
    
    Max.dist <- Mv.data %>% filter(dists == max(Mv.data$dists[Mv.data$case_ == TRUE]))
    
    Max.dist.step <- Max.dist$step_id_
    
    Mv.data <- Mv.data %>% mutate(inflec = ifelse(step_id_ <= Max.dist.step, 0, 1))
    
    all.data.inflec <- rbind(all.data.inflec, Mv.data)

  }
  
}

write.csv(all.data.inflec, "deerdata_4.csv")
