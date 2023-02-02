# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 08 - EHRM metrics (number per deer, distance, etc.)
# Author: Nathan D. Hooven
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Email: nathan.d.hooven@gmail.com
# Date started: 6 May 2022
# Date completed: 6 May 2022
# Date modified: 13 Jan 2023
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(sp)

#__________________________________________________________________________________________________
# 2. Read in data ----
#__________________________________________________________________________________________________

deer.data <- read.csv("deerdata_6.csv")

#__________________________________________________________________________________________________
# 3. Manipulate data for calculations ----
#__________________________________________________________________________________________________

deer.data.1 <- deer.data %>% filter(case_ == TRUE) %>%
                             dplyr::select(Deer, Sex, Age, MV.ID, burst_, x, y, t1_, step_id_, sl_, ta_, dists, inflec)

# number of excursions per deer
deer.n.exc <- deer.data.1 %>% group_by(Deer) %>%
                              summarize(n = length(unique(burst_)))

min(deer.n.exc$n)
max(deer.n.exc$n)
mean(deer.n.exc$n)
sd(deer.n.exc$n) / sqrt(nrow(deer.n.exc))

# furthest distance for excursions - min, max, median, mean
deer.dist.exc <- deer.data.1 %>% filter(grepl("EXC" , MV.ID)) %>%        # filter only excursions
                                 group_by(burst_) %>%                    # group by movement
                                 summarize(furthest.dist = max(dists))

min(deer.dist.exc$furthest.dist)
max(deer.dist.exc$furthest.dist)
median(deer.dist.exc$furthest.dist)
mean(deer.dist.exc$furthest.dist)
nrow(deer.dist.exc)

# furthest distance for dispersals - min, max, median, mean
deer.dist.disp <- deer.data.1 %>% filter(grepl("DISP" , MV.ID)) %>%        # filter only dispersals
                                  group_by(burst_) %>%                     # group by movement
                                  summarize(furthest.dist = max(dists))

min(deer.dist.disp$furthest.dist)
max(deer.dist.disp$furthest.dist)
median(deer.dist.disp$furthest.dist)
mean(deer.dist.disp$furthest.dist)
nrow(deer.dist.disp)

# plot these
deer.dist.exc$type <- "Excursions"
deer.dist.disp$type <- "Dispersals"

deer.dist.all <- rbind(deer.dist.exc, deer.dist.disp)

ggplot(deer.dist.all,
       aes(x = log(furthest.dist),
           color = type,
           fill = type,
           linetype = type)) +
  geom_density(alpha = 0.1,
               size = 1.05) +
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red")) +
  theme_bw() +
  xlab("ln(Furthest straight-line distance from first relocation) (m)") +
  ylab("Density") +
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.8, 0.8),
        panel.grid = element_blank())

# statistical test
shapiro.test(deer.dist.all$furthest.dist)

kruskal.test(formula = furthest.dist ~ type, data = deer.dist.all)

#__________________________________________________________________________________________________
# 4. Unfamiliar habitat as a function of distance ----
#__________________________________________________________________________________________________

# pivot
deer.data.pivot <- deer.data %>% pivot_longer(cols = 56:63) %>%
                                 filter(case_ == FALSE,
                                        name %in% c("D2.indiv.hr.175",
                                                    "D2.indiv.hr.250",
                                                    "D2.indiv.hr.350",
                                                    "D2.indiv.hr.500"))

deer.data.pivot$name <- factor(deer.data.pivot$name,
                               levels = c("D2.indiv.hr.175",
                                          "D2.indiv.hr.250",
                                          "D2.indiv.hr.350",
                                          "D2.indiv.hr.500"),
                               labels = c("175 m",
                                          "250 m",
                                          "350 m",
                                          "500 m"))

ggplot(data = deer.data.pivot, 
       aes(x = dists,
           y = value)) +
  facet_wrap(~name,
             scales = "free_y") +
  geom_point(alpha = 0.05) +
  geom_smooth(se = FALSE,
              color = "orange") +
  theme_bw() +
  xlab("Distance from starting relocation (m)") +
  ylab(expression(paste(D^2, " from indiv. HR locations"))) +
  theme(panel.grid = element_blank())
