# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 4 - Subsample pre-data for each EHRM
# Author: Nathan D. Hooven
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Email: nathan.hooven@uky.edu
# Date started: 9 Dec 2020
# Date completed: 9 Dec 2020
# Date modified: 5 May 2022
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(lubridate)     # work with dates
library(sp)            # work with spatial data

#__________________________________________________________________________________________________
# 2. Read in data ----
#__________________________________________________________________________________________________

all.data <- read.csv("deerdata_2.csv")

all.data$Timestamp <- as.POSIXct(ymd_hms(all.data$Timestamp, tz = "America/Chicago"))

EHRMs <- read.csv("EHRMs.csv")

EHRMs$Timestamp <- as.POSIXct(EHRMs$Timestamp, format = "%m/%d/%Y %H:%M", tz = "America/Chicago")

# remove all EHRM relocations from "pre HR" dataset
all.data.noEHRM <- all.data %>% anti_join(y = EHRMs, by = c("Deer", "UTME", "UTMN", "Timestamp"))

#__________________________________________________________________________________________________
# 3. Subset data for each EHRM for sampling ----
#__________________________________________________________________________________________________

all.pre.data <- data.frame()

for (i in unique(EHRMs$Deer)) {
  
  DeerID <- i
  
  Deer.data <- EHRMs %>% filter(Deer == DeerID)
  
  for (j in unique(Deer.data$MV.ID)) {
    
    MV <- j
    
    Mv.data <- Deer.data %>% filter(MV.ID == MV)
    
    # subset deer data for all (non-EHRM) locations before movement
    pre.data <- all.data.noEHRM %>% filter(Deer == DeerID & Timestamp < Mv.data$Timestamp[1])
    
    # add column with MV.ID
    pre.data$MV.ID <- MV
    
    # bind to master df
    all.pre.data <- rbind(all.pre.data, pre.data)
    
  }
  
}

#__________________________________________________________________________________________________
# 4. Write to .csv ----
#__________________________________________________________________________________________________

write.csv(all.pre.data, "all_pre_data.csv")
