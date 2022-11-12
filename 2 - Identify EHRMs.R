# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 2 - Identify EHRMs with Jacobsen et al.'s algorithm (highly modified)
# Author: Nathan D. Hooven
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Email: nathan.hooven@uky.edu
# Date started: 5 Dec 2020
# Date completed: 9 Dec 2020
# Date modified: 
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(sp)            # work with spatial data
library(adehabitatHR)  # create home ranges
library(adehabitatLT)  # create ltraj objects
library(rgdal)         # save shapefiles
library(raster)		     # 
library(maptools)	     #
library(rgeos)		     #
library(geosphere)	   # 
library(ggplot2)       # visualize dynamics in movement metrics
library(ggspatial)     # visualize maps

#__________________________________________________________________________________________________
# 2. Read in raw data ----
#__________________________________________________________________________________________________

deerdata <- read.csv("deerdata_2.csv")

deerdata$Timestamp <- as.POSIXct(deerdata$Timestamp)

#__________________________________________________________________________________________________
# 3. Identify EHRM for each animal with a for loop ----
#__________________________________________________________________________________________________

# create blank dataframe of relocations >= 500 m outside of the previous 30 day's HR isopleth
outside.relocations <- data.frame()

deer.names <- as.vector(unique(deerdata$Deer))

for (i in deer.names){
  
  DeerID <- i 
  
  Animal <- deerdata %>% filter(Deer == DeerID)
  
  # how many days worth of data do we have?
  total.days <- round(as.numeric(Animal[nrow(Animal), 9] - Animal[1, 9]))
  
  # define how many 30-day periods to look at
  periods <- (total.days - 31)
  
  # define and fill lists of periods (both PreHR and moving windows)
  PREStart <- Animal[1, 9] - (1*24*60*60)
  PREEnd <- PREStart + (29*24*60*60)
  PreS <- PREStart
  PreE <- PREEnd

  WinStart <- Animal[1, 9] + (29*24*60*60)
  WinEnd <- WinStart + (2*24*60*60)
  WinS <- WinStart
  WinE <- WinEnd
 
  for(i in 1:periods){
    PreS[i] <- PREStart + i*(24*60*60)
    PreE[i] <- PREEnd + i*(24*60*60)
  }
  
  for(i in 1:periods){
    WinS[i] <- WinStart + i*(24*60*60)
    WinE[i] <- WinEnd + i*(24*60*60)
  }
  
  # create blank lists and data frames to hold all data
  PrePd <- list()
  WinPd <- list()
  out.BBMM <- list()
  inwindow <- list()
  alloutside <- data.frame()
  
  # fit BBMM HRs for each PreHR, and determine distance of outside relocations in each moving window
  for (p in 1:periods){
    
    # define relocations for PreHR p
    PrePd[[p]] <- Animal %>% filter(Timestamp >= PreS[p] & Timestamp <= PreE[p])
    
    # only keep unique rows to prevent ltraj from throwing an error
    PrePd[[p]] <- PrePd[[p]] %>% distinct(Timestamp, .keep_all = TRUE)
    
    # fit trajectory to PreHR p relocations
    PrePd.ltraj <- as.ltraj(xy = PrePd[[p]][,c("UTME", "UTMN")],
                             date = PrePd[[p]]$Timestamp,
                             id = DeerID,
                             proj4string = CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
    
    # Sig1 using a conservative estimate of Sig2
    lik <- liker(PrePd.ltraj, sig2 = 20, rangesig1 = c(0, 50))
    
    # Estimate BBMM
    out.BBMM[[p]] <- kernelbb(PrePd.ltraj, 
                              sig1 = ifelse(substr(DeerID, 1, 1) == "1", lik$`1`$sig1, lik$`2`$sig1),
                              sig2 = 20, 
                              grid = 100,
                              extent = 1)
    
    # pull 95% contour, add CRS, and convert to longlat for later calculations
    contour <- getverticeshr(out.BBMM[[p]], percent = 95)
    contour@proj4string <- CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
    contour.longlat <- spTransform(contour, CRS("+proj=longlat"))

    # define relocations for moving window p
    WinPd[[p]] <- Animal %>% filter(Timestamp >= WinS[p] & Timestamp <= WinE[p])
    
    # promote moving window relocations to SPDF
    h <- as.data.frame(WinPd[p])
    
    h <- SpatialPointsDataFrame(coords = h[,c(7, 8)],
                                data = h,
                                proj4string = CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
    
    # create column of T/F values depending upon whether moving window relocations are within 95% BBMM
    h$inout <- !is.na(over(h, as(contour, "SpatialPolygons"))) 
    
    # create SP layer of UTM coordinates from moving window relocations, and convert to longlat
    utmcoor <- SpatialPoints(cbind(h$UTME,h$UTMN), proj4string = CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
    longlatcoor <- spTransform(utmcoor, CRS("+proj=longlat"))
    
    # calculate distances (in m) from moving window relocations and 95% BBMM isopleth
    h$dist <- dist2Line(longlatcoor, contour.longlat, distfun = distHaversine)
    
    k <- h@data
    k$period <- p
    inwindow[[p]] <- k
  }
  
  allinwindow <- do.call(rbind, inwindow)
  
  alloutside <- allinwindow[allinwindow$inout %in% FALSE,]
  
  alloutsidedist <- alloutside[alloutside$dist[,1] >= 500,]
  
  outside.relocations <- rbind(outside.relocations, alloutsidedist)
  
}

# write to csv
write.csv(outside.relocations, "outsiderelocations_final.csv")

# add ID column that actually makes sense
outside.relocations <- plyr::join(deerdata, outside.relocations, by = c("Deer", "UTME", "UTMN", "Timestamp"), type = "right")

names(outside.relocations)[11:15] <- c("X.2", "Sex.2", "Age.2", "Longitude.2", "Latitude.2")

outside.relocations <- outside.relocations %>% dplyr::select(X, Deer, Sex, Age, Longitude, Latitude, UTME, UTMN, Timestamp, inout, dist, period, ID)

#__________________________________________________________________________________________________
# 4. Confirm and delineate EHRM ----
#__________________________________________________________________________________________________

outside.relocations$Timestamp <- as.POSIXct(outside.relocations$Timestamp)

# define which animal and period we're looking at 
DeerID <- "172"
Period <- "207"

# subset the "outside.relocations" data frame to look at the movement itself
Movement <- outside.relocations %>% filter(Deer == DeerID & period == Period)

nrow(Movement)

# define the first relocation
First.reloc <- Movement[1,]

# subset the data before the movement began
Animal <- deerdata %>% filter(Deer == DeerID & Timestamp < First.reloc$Timestamp)

# create trajectory for movement
Movement.sp <- SpatialPoints(coords = Movement[,c("UTME", "UTMN")],
                             proj4string = CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

Movement <- Movement %>% distinct(Timestamp, .keep_all = TRUE)

Movement.ltraj <- as.ltraj(xy = Movement[,c("UTME", "UTMN")],
                           date = Movement$Timestamp,
                           id = DeerID,
                           proj4string = CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

Movement.traj <- ltraj2sldf(Movement.ltraj)

# create 95% BBMM HR for PreHR (all relocations, not just 30 days) 
# fit trajectory to PreHR relocations
Animal <- Animal %>% distinct(Timestamp, .keep_all = TRUE)

Animal.ltraj <- as.ltraj(xy = Animal[,c("UTME", "UTMN")],
                         date = Animal$Timestamp,
                         id = DeerID,
                         proj4string = CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

# Sig1 using a conservative estimate of Sig2
lik <- liker(Animal.ltraj, sig2 = 20, rangesig1 = c(0, 50))

# Estimate BBMM
Animal.BBMM <- kernelbb(Animal.ltraj, 
                        sig1 = ifelse(substr(DeerID, 1, 1) == "1", lik$`1`$sig1, lik$`2`$sig1),
                        sig2 = 20, 
                        grid = 100,
                        extent = 1)

# pull 95% contour, add CRS
Animal.contour <- getverticeshr(Animal.BBMM, percent = 95)
Animal.contour@proj4string <- CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

# ggplot
ggplot() +
  layer_spatial(Movement.traj) +
  layer_spatial(Movement.sp) +
  layer_spatial(Animal.contour) +
  annotation_spatial(Movement.traj, color = "#0066CC", size = 1) +
  annotation_spatial(Movement.sp, color = "#0066CC", size = 3) +
  annotation_spatial(Animal.contour, color = "darkgray", size = 1) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

#__________________________________________________________________________________________________
# 5. Subset and visualize final EHRMs ----
#__________________________________________________________________________________________________

# create blank dataframes to hold it all
#EHRMs <- data.frame()

DeerID <- "252"
StartID <- 200609 - 1
EndID <- 	200971 + 1

# Bound the complete EHRM
Final.EHRM <- deerdata %>% filter(Deer == DeerID & ID >= StartID & ID <= EndID)

# subset all the relocations before the movement
Final.PreHR <- deerdata %>% filter(Deer == DeerID & ID < StartID)

# create a trajectory for the EHRM
Final.EHRM.sp <- SpatialPoints(coords = Final.EHRM[,c("UTME", "UTMN")],
                               proj4string = CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

Final.EHRM <- Final.EHRM %>% distinct(Timestamp, .keep_all = TRUE)

Final.EHRM.ltraj <- as.ltraj(xy = Final.EHRM[,c("UTME", "UTMN")],
                             date = Final.EHRM$Timestamp,
                             id = DeerID,
                             proj4string = CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

Final.EHRM.traj <- ltraj2sldf(Final.EHRM.ltraj)

# create 95% BBMM HR for PreHR (all relocations, not just 30 days) 
# fit trajectory to PreHR relocations
Final.PreHR <- Final.PreHR %>% distinct(Timestamp, .keep_all = TRUE)

Final.PreHR.ltraj <- as.ltraj(xy = Final.PreHR[,c("UTME", "UTMN")],
                              date = Final.PreHR$Timestamp,
                              id = DeerID,
                              proj4string = CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

Final.PreHR.traj <- ltraj2sldf(Final.PreHR.ltraj)

# estimate 95% MCP
Final.MCP <- mcp(SpatialPoints(coords = Final.PreHR[,c("UTME", "UTMN")],
                               proj4string = CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")))

# ggplot
ggplot() +
  layer_spatial(Final.PreHR.traj) +
  layer_spatial(Final.EHRM.traj) +
  layer_spatial(Final.EHRM.sp) +
  annotation_spatial(Final.EHRM.traj, color = "#0066CC", size = 1) +
  annotation_spatial(Final.EHRM.sp, color = "black", size = 1.5) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

# create new MV.ID column and bind to master data frames
MV.ID <- "DISP"

# EHRMs
nrow(Final.EHRM)

Final.EHRM <- Final.EHRM %>% dplyr::mutate(MV.ID = MV.ID)

EHRMs <- rbind(EHRMs, Final.EHRM)

write.csv(EHRMs, "EHRMs.csv")

