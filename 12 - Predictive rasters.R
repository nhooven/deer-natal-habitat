# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 12 - Predictive rasters
# Author: Nathan D. Hooven
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Email: nathan.hooven@uky.edu
# Date started: 13 May 2022
# Date completed: 27 May 2022 
# Date modified: 
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(raster)        # work with rasters
library(sf)            # bounding box

#__________________________________________________________________________________________________
# 2. Load in models ----
#__________________________________________________________________________________________________

load("SSF_models.RData")

# pre-EHRM data
all.pre.data <- read.csv("all_pre_data_2.csv")

#__________________________________________________________________________________________________
# 3. Load in rasters ----
#__________________________________________________________________________________________________

# raster directory
raster.dir.dist <- "G:/Deer project/Deer project 5-2022/Rasters/Distance"

# forest distance
forest.dist <- raster(paste0(raster.dir.dist, "/forest_dist_clip.tif"))

# stream distance
stream.dist <- raster(paste0(raster.dir.dist, "/stream_dist.tif"))

# raster directory
raster.dir.frag <- "G:/Deer project/Deer project 5-2022/Rasters/Fragstats"

# 175-m buffer
pland.dev.175 <- raster(paste0(raster.dir.frag, "/175", "/pland_3.tif"))
pland.forest.175 <- raster(paste0(raster.dir.frag, "/175", "/pland_4.tif"))
pland.open.175 <- raster(paste0(raster.dir.frag, "/175", "/pland_5.tif"))
ai.175 <- raster(paste0(raster.dir.frag, "/175", "/ai.tif"))
iji.175 <- raster(paste0(raster.dir.frag, "/175", "/iji.tif"))
ed.175 <- raster(paste0(raster.dir.frag, "/175", "/ed.tif"))
prd.175 <- raster(paste0(raster.dir.frag, "/175", "/prd.tif"))

frag.175 <- stack(pland.dev.175, pland.forest.175, pland.open.175, 
                  ai.175, iji.175, ed.175, prd.175)

# 250-m buffer
pland.dev.250 <- raster(paste0(raster.dir.frag, "/250", "/pland_3.tif"))
pland.forest.250 <- raster(paste0(raster.dir.frag, "/250", "/pland_4.tif"))
pland.open.250 <- raster(paste0(raster.dir.frag, "/250", "/pland_5.tif"))
ai.250 <- raster(paste0(raster.dir.frag, "/250", "/ai.tif"))
iji.250 <- raster(paste0(raster.dir.frag, "/250", "/iji.tif"))
ed.250 <- raster(paste0(raster.dir.frag, "/250", "/ed.tif"))
prd.250 <- raster(paste0(raster.dir.frag, "/250", "/prd.tif"))

frag.250 <- stack(pland.dev.250, pland.forest.250, pland.open.250, 
                  ai.250, iji.250, ed.250, prd.250)

#__________________________________________________________________________________________________
# 4. Transform and scale distance rasters ----
#__________________________________________________________________________________________________
# 4a. Summer ----
#__________________________________________________________________________________________________

forest.dist.log <- log(forest.dist + 1)

forest.dist.summer <- (forest.dist.log - mean(deer.data.summer$dForest.P)) / sd(deer.data.summer$dForest.P)

stream.dist.summer <- (stream.dist - mean(deer.data.summer$dStream)) / sd(deer.data.summer$dStream)

#__________________________________________________________________________________________________
# 4b. Dispersal ----
#__________________________________________________________________________________________________

forest.dist.disp <- (forest.dist.log - mean(deer.data.disp$dForest.P)) / sd(deer.data.disp$dForest.P)

stream.dist.log <- log(stream.dist + 1)

stream.dist.disp <- (stream.dist.log - mean(deer.data.disp$dStream.P)) / sd(deer.data.disp$dStream.P)

#__________________________________________________________________________________________________
# 5. Mahalanobis distance from ALL ----
#__________________________________________________________________________________________________
# 5a. Summer - 250
#__________________________________________________________________________________________________

# determine dimensions of frag rasters (switch nrow and ncol)
ed.250
ncol.250 <- 2206
nrow.250 <- 2153

# all HR locations
all.hr.250 <- all.pre.data %>% dplyr::select(dev.250, forest.250, open.250, 
                                                 ai.250, iji.250, ed.250, prd.250)

# convert to matrix
frag.250.matrix <- as.matrix(frag.250)

# convert to df for replacing values
frag.250.matrix.df <- as.data.frame(frag.250.matrix)

# replace -999 in the first 19386 and final 17231 rows with NA
frag.250.matrix.df[1:19386, ] <- NA

frag.250.matrix.df[4732287:4749518, ] <- NA

# replace -999 with 0 (for pland) and 100 (for iji)
frag.250.matrix.df <- frag.250.matrix.df %>% mutate(pland_3 = replace(pland_3, pland_3 == -999, 0),
                                                    pland_4 = replace(pland_4, pland_4 == -999, 0),
                                                    pland_5 = replace(pland_5, pland_5 == -999, 0),
                                                    iji = replace(iji, iji == -999, 100))

# convert back to matrix
frag.250.matrix.1 <- as.matrix(frag.250.matrix.df)

# Mahalanobis distance
frag.250.D.ALL <- mahalanobis(frag.250.matrix.1, center = colMeans(all.hr.250), cov = cov(all.hr.250))

# convert to matrix
frag.250.D.ALL.matrix <- matrix(frag.250.D.ALL, nrow = nrow.250, ncol = ncol.250)

# replace all values > 300
frag.250.D.ALL.matrix[frag.250.D.ALL.matrix > 300] <- NA

# convert to raster
frag.250.D.raster <- raster(frag.250.D.ALL.matrix)

# plot to make sure it did what we intended
raster::plot(frag.250.D.raster)

# log-transform
frag.250.D.raster.log <- log(frag.250.D.raster)

# scale
frag.250.D.raster.summer <- (frag.250.D.raster.log - mean(deer.data.summer$D2.all.hr.250.P)) / sd(deer.data.summer$D2.all.hr.250.P)

# assign coordinates and a crs
frag.250.D.raster.summer@extent <- ed.250@extent
frag.250.D.raster.summer@crs <- ed.250@crs

#__________________________________________________________________________________________________
# 5b. Disp - 175
#__________________________________________________________________________________________________

# determine dimensions of frag rasters (switch nrow and ncol)
ed.175
ncol.175 <- 2206
nrow.175 <- 2153

# all HR locations
all.hr.175 <- all.pre.data %>% dplyr::select(dev.175, forest.175, open.175, 
                                             ai.175, iji.175, ed.175, prd.175)

# convert to matrix
frag.175.matrix <- as.matrix(frag.175)

# convert to df for replacing values
frag.175.matrix.df <- as.data.frame(frag.175.matrix)

# replace -999 in the first 19386 and final 17231 rows with NA
frag.175.matrix.df[1:19386, ] <- NA

frag.175.matrix.df[4732287:4749518, ] <- NA

# replace -999 with 0 (for pland) and 100 (for iji)
frag.175.matrix.df <- frag.175.matrix.df %>% mutate(pland_3 = replace(pland_3, pland_3 == -999, 0),
                                                    pland_4 = replace(pland_4, pland_4 == -999, 0),
                                                    pland_5 = replace(pland_5, pland_5 == -999, 0),
                                                    iji = replace(iji, iji == -999, 100))

# convert back to matrix
frag.175.matrix.1 <- as.matrix(frag.175.matrix.df)

# Mahalanobis distance
frag.175.D.ALL <- mahalanobis(frag.175.matrix.1, center = colMeans(all.hr.175), cov = cov(all.hr.175))

# convert to matrix
frag.175.D.ALL.matrix <- matrix(frag.175.D.ALL, nrow = nrow.175, ncol = ncol.175)

# replace all values > 300
frag.175.D.ALL.matrix[frag.175.D.ALL.matrix > 300] <- NA

# convert to raster
frag.175.D.raster <- raster(frag.175.D.ALL.matrix)

# plot to make sure it did what we intended
raster::plot(frag.175.D.raster)

# log-transform
frag.175.D.raster.log <- log(frag.175.D.raster)

# scale
frag.175.D.raster.disp <- (frag.175.D.raster.log - mean(deer.data.summer$D2.all.hr.175.P)) / sd(deer.data.summer$D2.all.hr.175.P)

# assign coordinates and a crs
frag.175.D.raster.disp@extent <- ed.175@extent
frag.175.D.raster.disp@crs <- ed.175@crs

#__________________________________________________________________________________________________
# 6. Predictive rasters (base + all) ----
#__________________________________________________________________________________________________
# 6a. Summer ----
#__________________________________________________________________________________________________

# make sure origins are the same
origin(forest.dist.summer) <- c(0, 0)
origin(stream.dist.summer) <- c(0, 0)
origin(frag.250.D.raster.summer) <- c(0, 0)

pred.base.summer <- exp(summer.models[[41]]$coefficients[1]*forest.dist.summer +
                        summer.models[[41]]$coefficients[2]*stream.dist.summer +
                        summer.models[[41]]$coefficients[3]*frag.250.D.raster.summer)

# fill in NAs
pred.base.summer.1 <- focal(x = pred.base.summer, 
                            w = matrix(1, nc = 3, nr = 3),
                            fun = mean, 
                            NAonly = T, 
                            na.rm = T)

# plot
raster::plot(pred.base.summer.1)

#__________________________________________________________________________________________________
# 6b. Disp ----
#__________________________________________________________________________________________________

# make sure origins are the same
origin(forest.dist.disp) <- c(0, 0)
origin(stream.dist.disp) <- c(0, 0)
origin(frag.175.D.raster.disp) <- c(0, 0)

# resample to reduce random streaking
frag.175.D.raster.disp.res <- resample(frag.175.D.raster.disp, stream.dist.disp)

pred.base.disp <- exp(disp.models[[41]]$coefficients[1]*forest.dist.disp +
                      disp.models[[41]]$coefficients[2]*stream.dist.disp +
                      disp.models[[41]]$coefficients[3]*frag.175.D.raster.disp.res)

# fill in NAs
pred.base.disp.1 <- focal(x = pred.base.disp, 
                          w = matrix(1, nc = 3, nr = 3),
                          fun = mean, 
                          NAonly = T, 
                          na.rm = T)

# plot
raster::plot(pred.base.disp.1)

#__________________________________________________________________________________________________
# 7. Mahalanobis distance from IND and predictive rasters ----
#__________________________________________________________________________________________________
# 7a. Summer - 175
#__________________________________________________________________________________________________

all.pre.data$burst_ <- paste0(all.pre.data$Deer, ".", all.pre.data$MV.ID)

# create a RasterStack to hold all rasters
summer.ind.stack <- stack()

for (x in unique(deer.data.summer.s$burst_)) {
  
  focal.mv <- x
  
  # indiv HR locations
  indiv.hr.175 <- all.pre.data %>% filter(burst_ == focal.mv) %>% 
                                   dplyr::select(dev.175, forest.175, open.175, 
                                                 ai.175, iji.175, ed.175, prd.175)
  
  # Mahalanobis distance
  frag.175.D.IND <- mahalanobis(frag.175.matrix.1, center = colMeans(indiv.hr.175), cov = cov(indiv.hr.175))
  
  # convert to matrix
  frag.175.D.IND.matrix <- matrix(frag.175.D.IND, nrow = nrow.175, ncol = ncol.175)
  
  # determine the D for the outer limits of the raster
  outer.limit <- frag.175.D.IND.matrix[1, 11]
  
  # replace all values of the outer limit with NA
  frag.175.D.IND.matrix[frag.175.D.IND.matrix == outer.limit] <- NA
  
  # convert to raster
  frag.175.D.IND.raster <- raster(frag.175.D.IND.matrix)
  
  # log-transform
  frag.175.D.IND.raster.log <- log(frag.175.D.IND.raster)
  
  # scale
  frag.175.D.IND.raster.summer <- (frag.175.D.IND.raster.log - mean(deer.data.summer$D2.indiv.hr.175.P)) / sd(deer.data.summer$D2.indiv.hr.175.P)
  
  # assign coordinates and a crs
  frag.175.D.IND.raster.summer@extent <- ed.175@extent
  frag.175.D.IND.raster.summer@crs <- ed.175@crs
  
  # make sure origins are the same
  origin(frag.175.D.IND.raster.summer) <- c(0, 0)
  
  pred.ind.summer <- exp(summer.models[[42]]$coefficients[1]*forest.dist.summer +
                         summer.models[[42]]$coefficients[2]*stream.dist.summer +
                         summer.models[[42]]$coefficients[3]*frag.250.D.raster.summer +
                         summer.models[[42]]$coefficients[4]*frag.175.D.IND.raster.summer)
  
  # fill in NAs
  pred.ind.summer.1 <- focal(x = pred.ind.summer, 
                              w = matrix(1, nc = 3, nr = 3),
                              fun = mean, 
                              NAonly = T, 
                              na.rm = T)
  
  # stack
  summer.ind.stack <- stack(summer.ind.stack, pred.ind.summer.1)
  
}

#__________________________________________________________________________________________________
# 7b. Disp - 175
#__________________________________________________________________________________________________

# create a RasterStack to hold all rasters
disp.ind.stack <- stack()

for (x in unique(deer.data.disp.s$burst_)) {
  
  focal.mv <- x
  
  # indiv HR locations
  indiv.hr.175 <- all.pre.data %>% filter(burst_ == focal.mv) %>% 
                                   dplyr::select(dev.175, forest.175, open.175, 
                                                 ai.175, iji.175, ed.175, prd.175)
  
  # Mahalanobis distance
  frag.175.D.IND <- mahalanobis(frag.175.matrix.1, center = colMeans(indiv.hr.175), cov = cov(indiv.hr.175))
  
  # convert to matrix
  frag.175.D.IND.matrix <- matrix(frag.175.D.IND, nrow = nrow.175, ncol = ncol.175)
  
  # determine the D for the outer limits of the raster
  outer.limit <- frag.175.D.IND.matrix[1, 11]
  
  # replace all values of the outer limit with NA
  frag.175.D.IND.matrix[frag.175.D.IND.matrix == outer.limit] <- NA
  
  # convert to raster
  frag.175.D.IND.raster <- raster(frag.175.D.IND.matrix)
  
  # log-transform
  frag.175.D.IND.raster.log <- log(frag.175.D.IND.raster)
  
  # scale
  frag.175.D.IND.raster.disp <- (frag.175.D.IND.raster.log - mean(deer.data.disp$D2.indiv.hr.175.P)) / sd(deer.data.disp$D2.indiv.hr.175.P)
  
  # assign coordinates and a crs
  frag.175.D.IND.raster.disp@extent <- ed.175@extent
  frag.175.D.IND.raster.disp@crs <- ed.175@crs
  
  # make sure origins are the same
  origin(frag.175.D.IND.raster.disp) <- c(0, 0)
  
  # resample to reduce random streaking
  frag.175.D.IND.raster.disp.res <- resample(frag.175.D.IND.raster.disp, stream.dist.disp)
  
  pred.ind.disp <- exp(disp.models[[42]]$coefficients[1]*forest.dist.disp +
                       disp.models[[42]]$coefficients[2]*stream.dist.disp +
                       disp.models[[42]]$coefficients[3]*frag.175.D.raster.disp +
                       disp.models[[42]]$coefficients[4]*frag.175.D.IND.raster.disp.res)
  
  # fill in NAs
  pred.ind.disp.1 <- focal(x = pred.ind.disp, 
                           w = matrix(1, nc = 3, nr = 3),
                           fun = mean, 
                           NAonly = T, 
                           na.rm = T)
  
  # stack
  disp.ind.stack <- stack(disp.ind.stack, pred.ind.disp.1)
  
}

#__________________________________________________________________________________________________
# 9. Discrete rasters ----

# extract bounding box
raster.bbox <- st_as_sfc(st_bbox(pred.base.summer.1))

#__________________________________________________________________________________________________
# 9a. Summer ----
#__________________________________________________________________________________________________

# define breaks
breaks.summer <- quantile(pred.base.summer.1, seq(0, 1, 1/10), na.rm = TRUE)

# cut into discrete categories
summer.base.discrete <- raster::cut(pred.base.summer.1, breaks = breaks.summer)

# convert to points
summer.base.discrete.pts <- rasterToPoints(summer.base.discrete, spatial = TRUE)

# convert to data frame
summer.base.discrete.df <- data.frame(summer.base.discrete.pts)

# convert to discrete data
summer.base.discrete.df$layer <- as.factor(summer.base.discrete.df$layer)

summer.base.plot <- ggplot() + 
                    geom_raster(data = summer.base.discrete.df, aes(x = x, 
                                                                    y = y, 
                                                                    fill = layer)) +
                    geom_sf(data = raster.bbox,
                            size = 1.5,
                            color = "black",
                            fill = NA) +
                    scale_fill_viridis_d(option = "plasma") +
                    theme_bw() +
                    theme(panel.grid = element_blank(),
                          legend.position = "none",
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          panel.border = element_blank(),
                          legend.background = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())

# example movements
# cut into discrete categories
summer.ind1.discrete <- raster::cut(summer.ind.stack$layer.1.1.1.1, breaks = breaks.summer)

# convert to points
summer.ind1.discrete.pts <- rasterToPoints(summer.ind1.discrete, spatial = TRUE)

# convert to data frame
summer.ind1.discrete.df <- data.frame(summer.ind1.discrete.pts)

# convert to discrete data
summer.ind1.discrete.df$layer <- as.factor(summer.ind1.discrete.df$layer)

summer.ind1.plot <- ggplot() + 
                    geom_raster(data = summer.base.discrete.df, aes(x = x, 
                                                                    y = y, 
                                                                    fill = layer)) +
                    geom_raster(data = summer.ind1.discrete.df, aes(x = x, 
                                                                    y = y, 
                                                                    fill = layer)) +
                    geom_sf(data = raster.bbox,
                            size = 1.5,
                            color = "black",
                            fill = NA) +
                    scale_fill_viridis_d(option = "plasma") +
                    theme_bw() +
                    theme(panel.grid = element_blank(),
                          legend.position = "none",
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          panel.border = element_blank(),
                          legend.background = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())

# calculate difference raster
summer.ind1.diff <- summer.ind1.discrete - summer.base.discrete

# convert to points
summer.ind1.diff.pts <- rasterToPoints(summer.ind1.diff, spatial = TRUE)

# convert to data frame
summer.ind1.diff.df <- data.frame(summer.ind1.diff.pts)

# convert to discrete data
summer.ind1.diff.df$layer <- as.factor(summer.ind1.diff.df$layer)

summer.ind1.diff.plot <- ggplot() + 
                    geom_raster(data = summer.ind1.diff.df, aes(x = x, 
                                                                    y = y, 
                                                                    fill = layer)) +
                    geom_sf(data = raster.bbox,
                            size = 1.5,
                            color = "black",
                            fill = NA) +
                    scale_fill_discrete(type = c(adjustcolor("purple4", alpha.f = 1.0), 
                                                 adjustcolor("purple4", alpha.f = 0.9), 
                                                 adjustcolor("purple4", alpha.f = 0.8), 
                                                 "white", 
                                                 adjustcolor("darkgreen", alpha.f = 0.6), 
                                                 adjustcolor("darkgreen", alpha.f = 0.7),
                                                 adjustcolor("darkgreen", alpha.f = 0.8),
                                                 adjustcolor("darkgreen", alpha.f = 0.9),
                                                 adjustcolor("darkgreen", alpha.f = 1.0))) +
                    theme_bw() +
                    theme(panel.grid = element_blank(),
                          legend.position = "none",
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          panel.border = element_blank(),
                          legend.background = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())

# plot together
cowplot::plot_grid(summer.base.plot, summer.ind1.plot, summer.ind1.diff.plot,
                   nrow = 1)

#__________________________________________________________________________________________________
# 9b. disp ----
#__________________________________________________________________________________________________

# define breaks
breaks.disp <- quantile(pred.base.disp.1, seq(0, 1, 1/10), na.rm = TRUE)

# cut into discrete categories
disp.base.discrete <- raster::cut(pred.base.disp.1, breaks = breaks.disp)

# convert to points
disp.base.discrete.pts <- rasterToPoints(disp.base.discrete, spatial = TRUE)

# convert to data frame
disp.base.discrete.df <- data.frame(disp.base.discrete.pts)

# convert to discrete data
disp.base.discrete.df$layer <- as.factor(disp.base.discrete.df$layer)

disp.base.plot <- ggplot() + 
                    geom_raster(data = disp.base.discrete.df, aes(x = x, 
                                                                    y = y, 
                                                                    fill = layer)) +
                    geom_sf(data = raster.bbox,
                            size = 1.5,
                            color = "black",
                            fill = NA) +
                    scale_fill_viridis_d(option = "plasma") +
                    theme_bw() +
                    theme(panel.grid = element_blank(),
                          legend.position = "none",
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          panel.border = element_blank(),
                          legend.background = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())

# example movements
# cut into discrete categories
disp.ind1.discrete <- raster::cut(disp.ind.stack$layer.2.2.1, breaks = breaks.disp)

# convert to points
disp.ind1.discrete.pts <- rasterToPoints(disp.ind1.discrete, spatial = TRUE)

# convert to data frame
disp.ind1.discrete.df <- data.frame(disp.ind1.discrete.pts)

# convert to discrete data
disp.ind1.discrete.df$layer <- as.factor(disp.ind1.discrete.df$layer)

disp.ind1.plot <- ggplot() + 
                    geom_raster(data = disp.base.discrete.df, aes(x = x, 
                                                                    y = y, 
                                                                    fill = layer)) +
                    geom_raster(data = disp.ind1.discrete.df, aes(x = x, 
                                                                    y = y, 
                                                                    fill = layer)) +
                    geom_sf(data = raster.bbox,
                            size = 1.5,
                            color = "black",
                            fill = NA) +
                    scale_fill_viridis_d(option = "plasma") +
                    theme_bw() +
                    theme(panel.grid = element_blank(),
                          legend.position = "none",
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          panel.border = element_blank(),
                          legend.background = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())

# calculate difference raster
disp.ind1.diff <- disp.ind1.discrete - disp.base.discrete

# convert to points
disp.ind1.diff.pts <- rasterToPoints(disp.ind1.diff, spatial = TRUE)

# convert to data frame
disp.ind1.diff.df <- data.frame(disp.ind1.diff.pts)

# convert to discrete data
disp.ind1.diff.df$layer <- as.factor(disp.ind1.diff.df$layer)

disp.ind1.diff.plot <- ggplot() + 
                    geom_raster(data = disp.ind1.diff.df, aes(x = x, 
                                                                    y = y, 
                                                                    fill = layer)) +
                    geom_sf(data = raster.bbox,
                            size = 1.5,
                            color = "black",
                            fill = NA) +
                    scale_fill_discrete(type = c(adjustcolor("purple4", alpha.f = 1.0), 
                                                 adjustcolor("purple4", alpha.f = 0.9), 
                                                 adjustcolor("purple4", alpha.f = 0.8), 
                                                 adjustcolor("purple4", alpha.f = 0.7), 
                                                 adjustcolor("purple4", alpha.f = 0.6), 
                                                 adjustcolor("purple4", alpha.f = 0.5), 
                                                 adjustcolor("purple4", alpha.f = 0.4), 
                                                 adjustcolor("purple4", alpha.f = 0.3),  
                                                 adjustcolor("purple4", alpha.f = 0.2), 
                                                 "white", 
                                                 adjustcolor("darkgreen", alpha.f = 0.3), 
                                                 adjustcolor("darkgreen", alpha.f = 0.4),
                                                 adjustcolor("darkgreen", alpha.f = 0.5),
                                                 adjustcolor("darkgreen", alpha.f = 0.6), 
                                                 adjustcolor("darkgreen", alpha.f = 0.7),
                                                 adjustcolor("darkgreen", alpha.f = 0.8),
                                                 adjustcolor("darkgreen", alpha.f = 0.9),
                                                 adjustcolor("darkgreen", alpha.f = 1.0))) +
                    theme_bw() +
                    theme(panel.grid = element_blank(),
                          legend.position = "none",
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          panel.border = element_blank(),
                          legend.background = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())

# plot together
cowplot::plot_grid(disp.base.plot, disp.ind1.plot, disp.ind1.diff.plot,
                   nrow = 1)

#__________________________________________________________________________________________________
# 9c. All plots ----
#__________________________________________________________________________________________________

cowplot::plot_grid(summer.base.plot, summer.ind1.plot, summer.ind1.diff.plot,
                   disp.base.plot, disp.ind1.plot, disp.ind1.diff.plot,
                   nrow = 2)

#__________________________________________________________________________________________________
# 10. Rasters for workflow diagram ----
#__________________________________________________________________________________________________

plot(forest.dist)

#__________________________________________________________________________________________________
# 11. Save image ----
#__________________________________________________________________________________________________

save.image("final_rasters.RData")
