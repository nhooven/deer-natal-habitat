# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 6 - Sample landscape variables
# Author: Nathan D. Hooven
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Email: nathan.hooven@uky.edu
# Date started: 5 May 2022
# Date completed: 5 May 2022
# Date modified: 7 Jun 2022
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(sp)            # work with spatial data
library(raster)        # sample from rasters

#__________________________________________________________________________________________________
# 2. Read in relocation datasets ----
#__________________________________________________________________________________________________

# pre-EHRM data
all.pre.data <- read.csv("all_pre_data.csv")

# EHRM data
ehrm.data <- read.csv("deerdata_4.csv")

# UTM projection
UTM.proj <- CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

#__________________________________________________________________________________________________
# 3. Read in rasters ----
#__________________________________________________________________________________________________
# 3a. Distance rasters ----
#__________________________________________________________________________________________________

# raster directory
raster.dir.dist <- "G:/Deer project/Deer project 5-2022/Rasters/Distance"

# ag distance
ag.dist <- raster(paste0(raster.dir.dist, "/ag_dist_clip.tif"))

# forest distance
forest.dist <- raster(paste0(raster.dir.dist, "/forest_dist_clip.tif"))

# road distance
road.dist <- raster(paste0(raster.dir.dist, "/road_dist_clip.tif"))

# stream distance
stream.dist <- raster(paste0(raster.dir.dist, "/stream_dist.tif"))

#__________________________________________________________________________________________________
# 3b. Fragstats rasters ----
#__________________________________________________________________________________________________

# raster directory
raster.dir.frag <- "G:/Deer project/Deer project 5-2022/Rasters/Fragstats"

# 175-m buffer
pland.dev.175 <- raster(paste0(raster.dir.frag, "/175", "/pland_3.tif"))
pland.forest.175 <- raster(paste0(raster.dir.frag, "/175", "/pland_4.tif"))
pland.open.175 <- raster(paste0(raster.dir.frag, "/175", "/pland_5.tif"))
pland.ag.175 <- raster(paste0(raster.dir.frag, "/175", "/pland_6.tif"))
ai.175 <- raster(paste0(raster.dir.frag, "/175", "/ai.tif"))
iji.175 <- raster(paste0(raster.dir.frag, "/175", "/iji.tif"))
ed.175 <- raster(paste0(raster.dir.frag, "/175", "/ed.tif"))
prd.175 <- raster(paste0(raster.dir.frag, "/175", "/prd.tif"))

frag.175 <- stack(pland.dev.175, pland.forest.175, pland.open.175, pland.ag.175,
                  ai.175, iji.175, ed.175, prd.175)

# 250-m buffer
pland.dev.250 <- raster(paste0(raster.dir.frag, "/250", "/pland_3.tif"))
pland.forest.250 <- raster(paste0(raster.dir.frag, "/250", "/pland_4.tif"))
pland.open.250 <- raster(paste0(raster.dir.frag, "/250", "/pland_5.tif"))
pland.ag.250 <- raster(paste0(raster.dir.frag, "/250", "/pland_6.tif"))
ai.250 <- raster(paste0(raster.dir.frag, "/250", "/ai.tif"))
iji.250 <- raster(paste0(raster.dir.frag, "/250", "/iji.tif"))
ed.250 <- raster(paste0(raster.dir.frag, "/250", "/ed.tif"))
prd.250 <- raster(paste0(raster.dir.frag, "/250", "/prd.tif"))

frag.250 <- stack(pland.dev.250, pland.forest.250, pland.open.250, pland.ag.250,
                  ai.250, iji.250, ed.250, prd.250)

# 350-m buffer
pland.dev.350 <- raster(paste0(raster.dir.frag, "/350", "/pland_3.tif"))
pland.forest.350 <- raster(paste0(raster.dir.frag, "/350", "/pland_4.tif"))
pland.open.350 <- raster(paste0(raster.dir.frag, "/350", "/pland_5.tif"))
pland.ag.350 <- raster(paste0(raster.dir.frag, "/350", "/pland_6.tif"))
ai.350 <- raster(paste0(raster.dir.frag, "/350", "/ai.tif"))
iji.350 <- raster(paste0(raster.dir.frag, "/350", "/iji.tif"))
ed.350 <- raster(paste0(raster.dir.frag, "/350", "/ed.tif"))
prd.350 <- raster(paste0(raster.dir.frag, "/350", "/prd.tif"))

frag.350 <- stack(pland.dev.350, pland.forest.350, pland.open.350, pland.ag.350,
                  ai.350, iji.350, ed.350, prd.350)

# 500-m buffer
pland.dev.500 <- raster(paste0(raster.dir.frag, "/500", "/pland_3.tif"))
pland.forest.500 <- raster(paste0(raster.dir.frag, "/500", "/pland_4.tif"))
pland.open.500 <- raster(paste0(raster.dir.frag, "/500", "/pland_5.tif"))
pland.ag.500 <- raster(paste0(raster.dir.frag, "/500", "/pland_6.tif"))
ai.500 <- raster(paste0(raster.dir.frag, "/500", "/ai.tif"))
iji.500 <- raster(paste0(raster.dir.frag, "/500", "/iji.tif"))
ed.500 <- raster(paste0(raster.dir.frag, "/500", "/ed.tif"))
prd.500 <- raster(paste0(raster.dir.frag, "/500", "/prd.tif"))

frag.500 <- stack(pland.dev.500, pland.forest.500, pland.open.500, pland.ag.500,
                  ai.500, iji.500, ed.500, prd.500)

#__________________________________________________________________________________________________
# 4. Extract values for all points ----
#__________________________________________________________________________________________________
# 4a. Pre-EHRM data ----
#__________________________________________________________________________________________________

# convert to SP
all.pre.data.SP <- SpatialPoints(coords = all.pre.data[ ,c("UTME", "UTMN")],
                                 proj4string = UTM.proj)

# distance rasters
all.pre.data$dAg <- extract(ag.dist, all.pre.data.SP)
all.pre.data$dForest <- extract(forest.dist, all.pre.data.SP)
all.pre.data$dRoad <- extract(road.dist, all.pre.data.SP)
all.pre.data$dStream <- extract(stream.dist, all.pre.data.SP)

# frag metrics
# 175-m
pre.frag.df.175 <- as.data.frame(extract(frag.175, all.pre.data.SP))

names(pre.frag.df.175) <- c("dev.175", "forest.175", "open.175", "ag.175",
                            "ai.175", "iji.175", "ed.175", "prd.175")

# replace -999 with 0 for pland variables, and 100 for iji
pre.frag.df.175 <- pre.frag.df.175 %>% mutate(dev.175 = replace(dev.175, dev.175 == -999, 0),
                                              forest.175 = replace(forest.175, forest.175 == -999, 0),
                                              open.175 = replace(open.175, open.175 == -999, 0), 
                                              ag.175 = replace(ag.175, ag.175 == -999, 0),
                                              iji.175 = replace(iji.175, iji.175 == -999, 100))

# bind to full dataset
all.pre.data <- cbind(all.pre.data, pre.frag.df.175)

# 250-m
pre.frag.df.250 <- as.data.frame(extract(frag.250, all.pre.data.SP))

names(pre.frag.df.250) <- c("dev.250", "forest.250", "open.250", "ag.250",
                            "ai.250", "iji.250", "ed.250", "prd.250")

# replace -999 with 0 for pland variables, and 100 for iji
pre.frag.df.250 <- pre.frag.df.250 %>% mutate(dev.250 = replace(dev.250, dev.250 == -999, 0),
                                              forest.250 = replace(forest.250, forest.250 == -999, 0),
                                              open.250 = replace(open.250, open.250 == -999, 0), 
                                              ag.250 = replace(ag.250, ag.250 == -999, 0),
                                              iji.250 = replace(iji.250, iji.250 == -999, 100))

# bind to full dataset
all.pre.data <- cbind(all.pre.data, pre.frag.df.250)

# 350-m
pre.frag.df.350 <- as.data.frame(extract(frag.350, all.pre.data.SP))

names(pre.frag.df.350) <- c("dev.350", "forest.350", "open.350", "ag.350",
                            "ai.350", "iji.350", "ed.350", "prd.350")

# replace -999 with 0 for pland variables, and 100 for iji
pre.frag.df.350 <- pre.frag.df.350 %>% mutate(dev.350 = replace(dev.350, dev.350 == -999, 0),
                                              forest.350 = replace(forest.350, forest.350 == -999, 0),
                                              open.350 = replace(open.350, open.350 == -999, 0), 
                                              ag.350 = replace(ag.350, ag.350 == -999, 0),
                                              iji.350 = replace(iji.350, iji.350 == -999, 100))

# bind to full dataset
all.pre.data <- cbind(all.pre.data, pre.frag.df.350)

# 500-m
pre.frag.df.500 <- as.data.frame(extract(frag.500, all.pre.data.SP))

names(pre.frag.df.500) <- c("dev.500", "forest.500", "open.500", "ag.500",
                            "ai.500", "iji.500", "ed.500", "prd.500")

# replace -999 with 0 for pland variables, and 100 for iji
pre.frag.df.500 <- pre.frag.df.500 %>% mutate(dev.500 = replace(dev.500, dev.500 == -999, 0),
                                              forest.500 = replace(forest.500, forest.500 == -999, 0),
                                              open.500 = replace(open.500, open.500 == -999, 0), 
                                              ag.500 = replace(ag.500, ag.500 == -999, 0),
                                              iji.500 = replace(iji.500, iji.500 == -999, 100))

# bind to full dataset
all.pre.data <- cbind(all.pre.data, pre.frag.df.500)

#__________________________________________________________________________________________________
# 4b. EHRM data ----
#__________________________________________________________________________________________________

# convert to SP
ehrm.data.SP <- SpatialPoints(coords = ehrm.data[ ,c("x", "y")],
                              proj4string = UTM.proj)

# distance rasters
ehrm.data$dAg <- extract(ag.dist, ehrm.data.SP)
ehrm.data$dForest <- extract(forest.dist, ehrm.data.SP)
ehrm.data$dRoad <- extract(road.dist, ehrm.data.SP)
ehrm.data$dStream <- extract(stream.dist, ehrm.data.SP)

# frag metrics
# 175-m
ehrm.frag.df.175 <- as.data.frame(extract(frag.175, ehrm.data.SP))

names(ehrm.frag.df.175) <- c("dev.175", "forest.175", "open.175", "ag.175",
                            "ai.175", "iji.175", "ed.175", "prd.175")

# replace -999 with 0 for pland variables, and 100 for iji
ehrm.frag.df.175 <- ehrm.frag.df.175 %>% mutate(dev.175 = replace(dev.175, dev.175 == -999, 0),
                                              forest.175 = replace(forest.175, forest.175 == -999, 0),
                                              open.175 = replace(open.175, open.175 == -999, 0), 
                                              ag.175 = replace(ag.175, ag.175 == -999, 0),
                                              iji.175 = replace(iji.175, iji.175 == -999, 100))

# bind to full dataset
ehrm.data <- cbind(ehrm.data, ehrm.frag.df.175)

# 250-m
ehrm.frag.df.250 <- as.data.frame(extract(frag.250, ehrm.data.SP))

names(ehrm.frag.df.250) <- c("dev.250", "forest.250", "open.250", "ag.250",
                            "ai.250", "iji.250", "ed.250", "prd.250")

# replace -999 with 0 for pland variables, and 100 for iji
ehrm.frag.df.250 <- ehrm.frag.df.250 %>% mutate(dev.250 = replace(dev.250, dev.250 == -999, 0),
                                              forest.250 = replace(forest.250, forest.250 == -999, 0),
                                              open.250 = replace(open.250, open.250 == -999, 0), 
                                              ag.250 = replace(ag.250, ag.250 == -999, 0),
                                              iji.250 = replace(iji.250, iji.250 == -999, 100))

# bind to full dataset
ehrm.data <- cbind(ehrm.data, ehrm.frag.df.250)

# 350-m
ehrm.frag.df.350 <- as.data.frame(extract(frag.350, ehrm.data.SP))

names(ehrm.frag.df.350) <- c("dev.350", "forest.350", "open.350", "ag.350",
                            "ai.350", "iji.350", "ed.350", "prd.350")

# replace -999 with 0 for pland variables, and 100 for iji
ehrm.frag.df.350 <- ehrm.frag.df.350 %>% mutate(dev.350 = replace(dev.350, dev.350 == -999, 0),
                                              forest.350 = replace(forest.350, forest.350 == -999, 0),
                                              open.350 = replace(open.350, open.350 == -999, 0), 
                                              ag.350 = replace(ag.350, ag.350 == -999, 0),
                                              iji.350 = replace(iji.350, iji.350 == -999, 100))

# bind to full dataset
ehrm.data <- cbind(ehrm.data, ehrm.frag.df.350)

# 500-m
ehrm.frag.df.500 <- as.data.frame(extract(frag.500, ehrm.data.SP))

names(ehrm.frag.df.500) <- c("dev.500", "forest.500", "open.500", "ag.500",
                            "ai.500", "iji.500", "ed.500", "prd.500")

# replace -999 with 0 for pland variables, and 100 for iji
ehrm.frag.df.500 <- ehrm.frag.df.500 %>% mutate(dev.500 = replace(dev.500, dev.500 == -999, 0),
                                              forest.500 = replace(forest.500, forest.500 == -999, 0),
                                              open.500 = replace(open.500, open.500 == -999, 0), 
                                              ag.500 = replace(ag.500, ag.500 == -999, 0),
                                              iji.500 = replace(iji.500, iji.500 == -999, 100))

# bind to full dataset
ehrm.data <- cbind(ehrm.data, ehrm.frag.df.500)

#__________________________________________________________________________________________________
# 5. Write to .csvs ----
#__________________________________________________________________________________________________

write.csv(all.pre.data, "all_pre_data_2.csv")

write.csv(ehrm.data, "deerdata_5.csv")

#__________________________________________________________________________________________________
# 6. Plot for supplementary ----
#__________________________________________________________________________________________________

# 175-m buffer
pland.dev.175 <- reclassify(pland.dev.175, cbind(-999, 0), right = FALSE)
pland.forest.175 <- reclassify(pland.forest.175, cbind(-999, 0), right = FALSE)
pland.open.175 <- reclassify(pland.open.175, cbind(-999, 0), right = FALSE)
iji.175 <- reclassify(iji.175, cbind(-999, NA), right = FALSE)
ai.175 <- reclassify(ai.175, cbind(-999, NA), right = FALSE)
ed.175 <- reclassify(ed.175, cbind(-999, NA), right = FALSE)
prd.175 <- reclassify(prd.175, cbind(-999, NA), right = FALSE)

frag.reclass.175 <- stack(pland.dev.175, pland.forest.175, pland.open.175,
                          ai.175, iji.175, ed.175, prd.175)

raster::plot(frag.reclass.175)

# 250-m buffer
pland.dev.250 <- reclassify(pland.dev.250, cbind(-999, 0), right = FALSE)
pland.forest.250 <- reclassify(pland.forest.250, cbind(-999, 0), right = FALSE)
pland.open.250 <- reclassify(pland.open.250, cbind(-999, 0), right = FALSE)
iji.250 <- reclassify(iji.250, cbind(-999, NA), right = FALSE)
ai.250 <- reclassify(ai.250, cbind(-999, NA), right = FALSE)
ed.250 <- reclassify(ed.250, cbind(-999, NA), right = FALSE)
prd.250 <- reclassify(prd.250, cbind(-999, NA), right = FALSE)

frag.reclass.250 <- stack(pland.dev.250, pland.forest.250, pland.open.250,
                          ai.250, iji.250, ed.250, prd.250)

raster::plot(frag.reclass.250)

# 350-m buffer
pland.dev.350 <- reclassify(pland.dev.350, cbind(-999, 0), right = FALSE)
pland.forest.350 <- reclassify(pland.forest.350, cbind(-999, 0), right = FALSE)
pland.open.350 <- reclassify(pland.open.350, cbind(-999, 0), right = FALSE)
iji.350 <- reclassify(iji.350, cbind(-999, NA), right = FALSE)
ai.350 <- reclassify(ai.350, cbind(-999, NA), right = FALSE)
ed.350 <- reclassify(ed.350, cbind(-999, NA), right = FALSE)
prd.350 <- reclassify(prd.350, cbind(-999, NA), right = FALSE)

frag.reclass.350 <- stack(pland.dev.350, pland.forest.350, pland.open.350,
                          ai.350, iji.350, ed.350, prd.350)

raster::plot(frag.reclass.350)

# 500-m buffer
pland.dev.500 <- reclassify(pland.dev.500, cbind(-999, 0), right = FALSE)
pland.forest.500 <- reclassify(pland.forest.500, cbind(-999, 0), right = FALSE)
pland.open.500 <- reclassify(pland.open.500, cbind(-999, 0), right = FALSE)
iji.500 <- reclassify(iji.500, cbind(-999, NA), right = FALSE)
ai.500 <- reclassify(ai.500, cbind(-999, NA), right = FALSE)
ed.500 <- reclassify(ed.500, cbind(-999, NA), right = FALSE)
prd.500 <- reclassify(prd.500, cbind(-999, NA), right = FALSE)

frag.reclass.500 <- stack(pland.dev.500, pland.forest.500, pland.open.500,
                          ai.500, iji.500, ed.500, prd.500)

raster::plot(frag.reclass.500)

# base covariates
par(mfrow = c(1, 4))

raster::plot(ag.dist, main = "ag.dist")
raster::plot(forest.dist, main = "forest.dist")
raster::plot(stream.dist, main = "stream.dist")
raster::plot(road.dist, main = "road.dist")