# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 03 - Generate random steps and summarize movements
# Author: Nathan D. Hooven
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Email: nathan.d.hooven@gmail.com
# Date started: 9 Dec 2020
# Date completed: 9 Dec 2020
# Date modified: 7 Jun 2022
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(sp)            # work with spatial data
library(amt)           # generate random steps
library(circular)      # von mises distribution
library(scales)        # radians scale

#__________________________________________________________________________________________________
# 2. Read in raw data ----
#__________________________________________________________________________________________________

EHRMs <- read.csv("EHRMs.csv")

EHRMs$Timestamp <- as.POSIXct(EHRMs$Timestamp, format = "%m/%d/%Y %H:%M", tz = "America/Chicago")

#__________________________________________________________________________________________________
# 3. Make tracks for all data and create steps by burst ----
#__________________________________________________________________________________________________

all.tracks <- data.frame()

for (i in unique(EHRMs$Deer)) {
  
  DeerID <- i
  
  Deer.data <- EHRMs %>% filter(Deer == DeerID)
  
  for (j in unique(Deer.data$MV.ID)) {
    
    MV <- j
    
    Mv.data <- Deer.data %>% filter(MV.ID == MV)
    
    deer.track <- make_track(Mv.data,
                             .x = UTME,
                             .y = UTMN,
                             .t = Timestamp,
                             all_cols = TRUE,
                             check_duplicates = TRUE,
                             crs = sp::CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
    
    all.tracks <- rbind(all.tracks, deer.track)
    
  }
  
}

# create "burst" column
all.tracks <- all.tracks %>% mutate(burst_ = paste0(Deer, ".", MV.ID))

# create steps by burst
all.steps <- steps_by_burst(all.tracks, keep_cols = "start")

#__________________________________________________________________________________________________
# 4. Summarize step lengths ----
#__________________________________________________________________________________________________

mean(all.steps$sl_)

median(all.steps$sl_)

max(all.steps$sl_)

ggplot(data = all.steps, aes(x = sl_)) +
  geom_density(fill = "lightgray") +
  theme_bw() +
  geom_vline(xintercept = mean(all.steps$sl_, na.rm = TRUE), color = "blue") +
  geom_vline(xintercept = median(all.steps$sl_, na.rm = TRUE), color = "red") +
  geom_vline(xintercept = max(all.steps$sl_, na.rm = TRUE), color = "darkgreen")

#__________________________________________________________________________________________________
# 5. Generate random steps ----
#__________________________________________________________________________________________________

# determine optimal distributions for step lengths
# step lengths
dist.gam <- fit_distr(all.steps$sl_, "gamma")
dist.exp <- fit_distr(all.steps$sl_, "exp")

ggplot(data = all.steps, aes(x = sl_)) +
       theme_bw() +
       geom_density(size = 1.5) +
       stat_function(fun = dgamma, 
                     args = list(shape = dist.gam$params$shape, 
                                 scale = dist.gam$params$scale),
                     color = "red",
                     size = 1.5) +
       stat_function(fun = dexp, 
                     args = list(rate = dist.exp$params$rate),
                     color = "blue",
                     size = 1.5)

dist.vm <- fit_distr(all.steps$ta_, "vonmises")

# sample 50x the used relocations
random.steps <- random_steps(all.steps, 
                             n_control = 50,
                             rand_sl = random_numbers(dist.gam, n = 1e+05),
                             rand_ta = random_numbers(dist.vm, n = 1e+05))

# keep x1_ and x2_ for case_ == TRUE and keep x2_ and y2_ for case == FALSE
correct.xy <- random.steps %>% mutate(x = ifelse(case_ == TRUE, x1_, x2_),
                                      y = ifelse(case_ == TRUE, y1_, y2_))

# keep only the columns we need
final.steps <- correct.xy %>% select(Deer, Sex, Age, ID, MV.ID, burst_, x, y, t1_, step_id_, case_, sl_, ta_)

# write to csv
write.csv(final.steps, "deerdata_3.csv")

#__________________________________________________________________________________________________
# 6. Summarize each movement (number of relocations, number of steps, time elapsed) ----
#__________________________________________________________________________________________________

burst.summaries <- all.tracks %>% dplyr::group_by(burst_) %>%
                                  dplyr::summarize(relocations = n(),
                                                   steps = n() - 1)
Time.summary <- data.frame()

# create time elapsed variable
for (i in unique(EHRMs$Deer)) {
  
  DeerID <- i
  
  Deer.data <- EHRMs %>% filter(Deer == DeerID)
  
  for (j in unique(Deer.data$MV.ID)) {
    
    MV <- j
    
    Mv.data <- Deer.data %>% filter(MV.ID == MV)
    
    Deer.summary <- data.frame(burst_ = paste0(DeerID, ".", MV), 
                               Time = difftime(max(Mv.data$Timestamp), min(Mv.data$Timestamp), units = "hours"))
    
    Time.summary <- rbind(Time.summary, Deer.summary)
    
  }
  
}

burst.summaries <- plyr::join(burst.summaries, Time.summary, by = "burst_", type = "right")

write.csv(burst.summaries, "burst_summaries.csv")

#__________________________________________________________________________________________________
# 7. Plot distributions ----
#__________________________________________________________________________________________________

# gamma for step lengths
dist.gam.pts <- rgamma(100000, shape = dist.gam$params$shape, scale = dist.gam$params$scale)

dist.gam.den <- density(dist.gam.pts)

sl.plot <- ggplot() +
  theme_bw() +
  geom_line(aes(x = dist.gam.den$x,
                y = dist.gam.den$y),
            size = 1,
            color = "darkblue") +
  xlab("Step length (m)") +
  ylab("") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  coord_cartesian(xlim = c(0, quantile(all.steps$sl_, probs = 0.99)))

# Von Mises for turn angles
dist.vm.pts <- rvonmises(10, 
                         mu = dist.vm$params$mu, 
                         kappa = dist.vm$params$kappa, 
                         control.circular = list(units = "radians"))

dist.vm.den <- density.circular(all.steps$ta_,
                                kernel = c("vonmises"),
                                na.rm = TRUE,
                                bw = dist.vm$params$kappa,
                                control.circular = list(units = "radians"))

# create radians scale
pi_scales <- math_format(.x * pi, format = function(x) x / pi)

ta.plot <- ggplot() +
  theme_bw() +
  geom_line(aes(x = dist.vm.den$x,
                y = dist.vm.den$y),
            size = 1,
            color = "red") +
  xlab("Turn angle (rad)") +
  ylab("") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  scale_x_continuous(breaks = c(0, pi/2, pi, (3*pi)/2, 2*pi),
                     labels = pi_scales)

# plot together
cowplot::plot_grid(sl.plot, ta.plot,
                   rel_widths = c(0.9, 1))
            