# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 9 - Prediction plots
# Author: Nathan D. Hooven
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Email: nathan.hooven@uky.edu
# Date started: 6 May 2022
# Date completed: 6 May 2022
# Date modified: 7 Jun 2022
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(lubridate)     # work with dates
library(survival)      # conditional logistic regressions
library(AICcmodavg)    # model selection
library(amt)           # calculate RSS
library(broom)
library(cowplot)       # multiple plots
library(glmmTMB)       # mixed models

#__________________________________________________________________________________________________
# 2. Write parameter tables ----
#__________________________________________________________________________________________________

load("SSF_models.RData")

all.mod.tidy <- data.frame() 
  
# loop through models
for (x in 40:42) {
  
  which.mod <- x

  # spring
  mod.tidy.spring <- tidy(spring.models[[which.mod]]) %>% 
                     mutate(model = which.mod,
                            season = "spring")
  
  # summer
  mod.tidy.summer <- tidy(summer.models[[which.mod]]) %>% 
                     mutate(model = which.mod,
                            season = "summer")
  
  # fall
  mod.tidy.fall <- tidy(fall.models[[which.mod]]) %>% 
                     mutate(model = which.mod,
                            season = "fall")
  
  # disp
  mod.tidy.disp <- tidy(disp.models[[which.mod]]) %>% 
                     mutate(model = which.mod,
                            season = "disp")
  
  all.mod.tidy <- rbind(all.mod.tidy, mod.tidy.spring, mod.tidy.summer,
                        mod.tidy.fall, mod.tidy.disp)
    
}

# round estimates and SEs, and paste together
all.mod.tidy$estimate.rd <- round(all.mod.tidy$estimate, digits = 2)
all.mod.tidy$se.rd <- round(all.mod.tidy$std.error, digits = 2)

all.mod.tidy$est.se <- paste0(all.mod.tidy$estimate.rd, " (", all.mod.tidy$se.rd, ")")

# select columns and pivot
all.mod.sum <- all.mod.tidy %>% dplyr::select(term, model, est.se, season) %>% 
                                pivot_wider(names_from = "term",
                                            values_from = "est.se")

write.table(all.mod.sum, "clipboard", sep = "\t")

#__________________________________________________________________________________________________
# 3. Relative selection strength prediction plots ----
#__________________________________________________________________________________________________
# 3a. Spring ----
#__________________________________________________________________________________________________

# refit models as "fit_clogit"
spring.mod1 <- fit_clogit(case_ ~
                               dForest.P +
                               dStream + dStream.Q +
                               strata(step_id_),
                             data = deer.data.spring.s,
                          model = TRUE)

spring.mod2 <- fit_clogit(case_ ~
                               dForest.P +
                               dStream + dStream.Q +
                                D2.all.hr.175.P +
                               strata(step_id_),
                             data = deer.data.spring.s,
                          model = TRUE)

spring.mod3 <- fit_clogit(case_ ~
                               dForest.P +
                               dStream + dStream.Q +
                                D2.all.hr.175.P +
                               D2.indiv.hr.175.P +
                               strata(step_id_),
                             data = deer.data.spring.s,
                          model = TRUE)

spring.x2 <- data.frame(dForest.P = 0,
                        dStream = 0,
                        dStream.Q = 0,
                        D2.all.hr.175.P = 0,
                        D2.indiv.hr.175.P = 0)

spring.x1 <- data.frame(dForest.P = 0,
                        dStream = 0,
                        dStream.Q = 0,
                        D2.all.hr.175.P = 0,
                        D2.indiv.hr.175.P = seq(min(deer.data.spring.s$D2.indiv.hr.175.P),
                                              max(deer.data.spring.s$D2.indiv.hr.175.P),
                                              length.out = 100))

# base model 
spring.mod1.rss <- log_rss(spring.mod1,
                           x1 = spring.x1,
                           x2 = spring.x2,
                           ci = "se",
                           ci_level = 0.90)

# base + all model 
spring.mod2.rss <- log_rss(spring.mod2,
                           x1 = spring.x1,
                           x2 = spring.x2,
                           ci = "se",
                           ci_level = 0.90)

# base + all + ind model 
spring.mod3.rss <- log_rss(spring.mod3,
                           x1 = spring.x1,
                           x2 = spring.x2,
                           ci = "se",
                           ci_level = 0.90)

# calculate model-averaged predictions
spring.rss <- data.frame(rss_mod1 = spring.mod1.rss$df$log_rss*0.56,
                         lwr_mod1 = spring.mod1.rss$df$lwr*0.56,
                         upr_mod1 = spring.mod1.rss$df$upr*0.56,
                         rss_mod2 = spring.mod2.rss$df$log_rss*0.40,
                         lwr_mod2 = spring.mod2.rss$df$lwr*0.40,
                         upr_mod2 = spring.mod2.rss$df$upr*0.40,
                         rss_mod3 = spring.mod3.rss$df$log_rss*0.04,
                         lwr_mod3 = spring.mod3.rss$df$lwr*0.04,
                         upr_mod3 = spring.mod3.rss$df$upr*0.04)

spring.rss.sum <- data.frame(rss = spring.rss$rss_mod1 + spring.rss$rss_mod2 + spring.rss$rss_mod3,
                             lwr = spring.rss$lwr_mod1 + spring.rss$lwr_mod2 + spring.rss$lwr_mod3,
                             upr = spring.rss$upr_mod1 + spring.rss$upr_mod2 + spring.rss$upr_mod3,
                             x.range = spring.mod1.rss$df$D2.indiv.hr.175.P_x1)

# determine range for back-transformation
spring.scaled.range <- seq(min(deer.data.spring.s$D2.indiv.hr.175.P),
                           max(deer.data.spring.s$D2.indiv.hr.175.P),
                           length.out = 100)

spring.rss.sum$D2.indiv.hr.175_x1_back <- (spring.scaled.range*sd(deer.data.spring$D2.indiv.hr.175.P)) + mean(deer.data.spring$D2.indiv.hr.175.P)

# plot prediction
spring.rss.plot <- ggplot(data = spring.rss.sum, aes(x = exp(D2.indiv.hr.175_x1_back), y = rss)) +
       theme_bw() +
       geom_hline(yintercept = 0) +
       geom_vline(xintercept = exp(mean(deer.data.spring$D2.indiv.hr.175.P)), linetype = "dashed") +
       geom_line(size = 1.5, color = "darkgreen") +
       geom_ribbon(aes(ymin = lwr,
                       ymax = upr),
                   fill = "darkgreen",
                   alpha = 0.1) +
       theme(panel.grid = element_blank()) +
       ylab("ln(Relative selection strength)") +
       xlab(expression(paste(D^2, " from indiv. HR locations (175-m)"))) +
       coord_cartesian(xlim = c(0, quantile(deer.data.spring$D2.indiv.hr.175, 0.99)),
                       ylim = c(-0.01, 0.04)) +
       ggtitle("a")

#__________________________________________________________________________________________________
# 3b. Summer ----
#__________________________________________________________________________________________________

# refit models as "fit_clogit"
summer.mod1 <- fit_clogit(case_ ~
                               dForest.P +
                               dStream +
                               strata(step_id_),
                             data = deer.data.summer.s,
                          model = TRUE)

summer.mod2 <- fit_clogit(case_ ~
                               dForest.P +
                               dStream +
                               D2.all.hr.250.P +
                               strata(step_id_),
                             data = deer.data.summer.s,
                          model = TRUE)

summer.mod3 <- fit_clogit(case_ ~
                               dForest.P +
                               dStream +
                               D2.all.hr.250.P +
                               D2.indiv.hr.175.P +
                               strata(step_id_),
                             data = deer.data.summer.s,
                          model = TRUE)

summer.x1 <- data.frame(dForest.P = 0,
                        dStream = 0,
                        dStream.Q = 0,
                        D2.all.hr.250.P = 0,
                        D2.indiv.hr.175.P = seq(min(deer.data.summer.s$D2.indiv.hr.175.P),
                                              max(deer.data.summer.s$D2.indiv.hr.175.P),
                                              length.out = 100))

summer.x2 <- data.frame(dForest.P = 0,
                        dStream = 0,
                        dStream.Q = 0,
                        D2.all.hr.250.P = 0,
                        D2.indiv.hr.175.P = 0)

# base
summer.mod1.rss <- log_rss(summer.mod1,
                           x1 = summer.x1,
                           x2 = summer.x2,
                           ci = "se",
                           ci_level = 0.90)

# base + all
summer.mod2.rss <- log_rss(summer.mod2,
                           x1 = summer.x1,
                           x2 = summer.x2,
                           ci = "se",
                           ci_level = 0.90)

# base + all + ind model
summer.mod3.rss <- log_rss(summer.mod3,
                           x1 = summer.x1,
                           x2 = summer.x2,
                           ci = "se",
                           ci_level = 0.90)

# calculate model-averaged predictions
summer.rss <- data.frame(rss_mod1 = summer.mod1.rss$df$log_rss*0.45,
                         lwr_mod1 = summer.mod1.rss$df$lwr*0.45,
                         upr_mod1 = summer.mod1.rss$df$upr*0.45,
                         rss_mod2 = summer.mod2.rss$df$log_rss*0.20,
                         lwr_mod2 = summer.mod2.rss$df$lwr*0.20,
                         upr_mod2 = summer.mod2.rss$df$upr*0.20,
                         rss_mod3 = summer.mod3.rss$df$log_rss*0.35,
                         lwr_mod3 = summer.mod3.rss$df$lwr*0.35,
                         upr_mod3 = summer.mod3.rss$df$upr*0.35)

summer.rss.sum <- data.frame(rss = summer.rss$rss_mod1 + summer.rss$rss_mod2 + summer.rss$rss_mod3,
                             lwr = summer.rss$lwr_mod1 + summer.rss$lwr_mod2 + summer.rss$lwr_mod3,
                             upr = summer.rss$upr_mod1 + summer.rss$upr_mod2 + summer.rss$upr_mod3,
                             x.range = summer.mod1.rss$df$D2.indiv.hr.175.P_x1)

# determine range for back-transformation
summer.scaled.range <- seq(min(deer.data.summer.s$D2.indiv.hr.175.P),
                                              max(deer.data.summer.s$D2.indiv.hr.175.P),
                                              length.out = 100)

summer.rss.sum$x.range.back <- (summer.scaled.range*sd(deer.data.summer$D2.indiv.hr.175.P)) + mean(deer.data.summer$D2.indiv.hr.175.P)

# plot prediction
summer.rss.plot <- ggplot(data = summer.rss.sum, aes(x = exp(x.range.back), y = rss)) +
       theme_bw() +
       geom_hline(yintercept = 0) +
       geom_vline(xintercept = exp(mean(deer.data.summer$D2.indiv.hr.175.P)), linetype = "dashed") +
       geom_line(size = 1.5, color = "orange") +
       geom_ribbon(aes(ymin = lwr,
                       ymax = upr),
                   fill = "orange",
                   alpha = 0.1) +
       theme(panel.grid = element_blank()) +
       ylab("") +
       xlab(expression(paste(D^2, " from indiv. HR locations (175-m)"))) +
       coord_cartesian(xlim = c(0, quantile(deer.data.summer$D2.indiv.hr.175, 0.99)),
                       ylim = c(-0.35, 0.15)) +
       ggtitle("b")

#__________________________________________________________________________________________________
# 3c. Fall ----
#__________________________________________________________________________________________________

# refit models as "fit_clogit"
fall.mod1 <- fit_clogit(case_ ~
                               dAg.P +
                                dForest.P +
                                dRoad.P +
                               dStream +
                               strata(step_id_),
                             data = deer.data.fall.s,
                          model = TRUE)

fall.mod2 <- fit_clogit(case_ ~
                               dAg.P +
                                dForest.P +
                                dRoad.P +
                               dStream +
                               D2.all.hr.175.P +
                               strata(step_id_),
                             data = deer.data.fall.s,
                          model = TRUE)

fall.mod3 <- fit_clogit(case_ ~
                               dAg.P +
                                dForest.P +
                                dRoad.P +
                               dStream +
                               D2.all.hr.175.P +
                               D2.indiv.hr.350.P +
                               strata(step_id_),
                             data = deer.data.fall.s,
                          model = TRUE)

fall.x1 <- data.frame(dAg.P = 0,
                       dForest.P = 0, 
                       dRoad.P = 0,
                      dStream = 0,
                        D2.all.hr.175.P = 0,
                        D2.indiv.hr.350.P = seq(min(deer.data.fall.s$D2.indiv.hr.350.P),
                                              max(deer.data.fall.s$D2.indiv.hr.350.P),
                                              length.out = 100))

fall.x2 <- data.frame(dAg.P = 0,
                       dForest.P = 0, 
                       dRoad.P = 0,
                      dStream = 0,
                        D2.all.hr.175.P = 0,
                        D2.indiv.hr.350.P = 0)

# base
fall.mod1.rss <- log_rss(fall.mod1,
                           x1 = fall.x1,
                           x2 = fall.x2,
                           ci = "se",
                           ci_level = 0.90)

# base + all
fall.mod2.rss <- log_rss(fall.mod2,
                           x1 = fall.x1,
                           x2 = fall.x2,
                           ci = "se",
                           ci_level = 0.90)

# base + all + ind model
fall.mod3.rss <- log_rss(fall.mod3,
                           x1 = fall.x1,
                           x2 = fall.x2,
                           ci = "se",
                           ci_level = 0.90)

# calculate model-averaged predictions
fall.rss <- data.frame(rss_mod1 = fall.mod1.rss$df$log_rss*0.45,
                         lwr_mod1 = fall.mod1.rss$df$lwr*0.45,
                         upr_mod1 = fall.mod1.rss$df$upr*0.45,
                         rss_mod2 = fall.mod2.rss$df$log_rss*0.32,
                         lwr_mod2 = fall.mod2.rss$df$lwr*0.32,
                         upr_mod2 = fall.mod2.rss$df$upr*0.32,
                         rss_mod3 = fall.mod3.rss$df$log_rss*0.23,
                         lwr_mod3 = fall.mod3.rss$df$lwr*0.23,
                         upr_mod3 = fall.mod3.rss$df$upr*0.23)

fall.rss.sum <- data.frame(rss = fall.rss$rss_mod1 + fall.rss$rss_mod2 + fall.rss$rss_mod3,
                             lwr = fall.rss$lwr_mod1 + fall.rss$lwr_mod2 + fall.rss$lwr_mod3,
                             upr = fall.rss$upr_mod1 + fall.rss$upr_mod2 + fall.rss$upr_mod3,
                             x.range = fall.mod1.rss$df$D2.indiv.hr.350.P_x1)

# determine range for back-transformation
fall.scaled.range <- seq(min(deer.data.fall.s$D2.indiv.hr.350.P),
                                              max(deer.data.fall.s$D2.indiv.hr.350.P),
                                              length.out = 100)

fall.rss.sum$x.range.back <- (fall.scaled.range*sd(deer.data.fall$D2.indiv.hr.350.P)) + mean(deer.data.fall$D2.indiv.hr.350.P)

# plot prediction
fall.rss.plot <- ggplot(data = fall.rss.sum, aes(x = exp(x.range.back), y = rss)) +
       theme_bw() +
       geom_hline(yintercept = 0) +
       geom_vline(xintercept = 35, linetype = "dashed") +
       geom_line(size = 1.5, color = "red") +
       geom_ribbon(aes(ymin = lwr,
                       ymax = upr),
                   fill = "red",
                   alpha = 0.1) +
       theme(panel.grid = element_blank()) +
       ylab("") +
       xlab(expression(paste(D^2, " from indiv. HR locations (350-m)"))) +
       coord_cartesian(xlim = c(0, quantile(deer.data.fall$D2.indiv.hr.350, 0.99)),
                       ylim = c(-0.15, 0.15)) +
       scale_y_continuous(breaks = seq(-0.20, 0.15, 0.05)) +
       ggtitle("c")

#__________________________________________________________________________________________________
# 3d. Dispersal ----
#__________________________________________________________________________________________________

# refit model as "fit_clogit"
disp.mod1 <- fit_clogit(case_ ~
                               dForest.P +
                               dStream.P +
                               D2.all.hr.175.P +
                               D2.indiv.hr.175.P +
                               strata(step_id_),
                             data = deer.data.disp.s,
                          model = TRUE)

disp.x1 <- data.frame(dForest.P = 0,
                        dStream.P = 0,
                        D2.all.hr.175.P = 0,
                        D2.indiv.hr.175.P = seq(min(deer.data.disp.s$D2.indiv.hr.175.P),
                                              max(deer.data.disp.s$D2.indiv.hr.175.P),
                                              length.out = 100))

disp.x2 <- data.frame(dForest.P = 0,
                        dStream.P = 0,
                        D2.all.hr.175.P = 0,
                        D2.indiv.hr.175.P = 0)

# base + all + ind model
disp.mod1.rss <- log_rss(disp.mod1,
                           x1 = disp.x1,
                           x2 = disp.x2,
                           ci = "se",
                           ci_level = 0.90)

# determine range for back-transformation
disp.scaled.range <- seq(min(deer.data.disp.s$D2.indiv.hr.175.P),
                         max(deer.data.disp.s$D2.indiv.hr.175.P),
                         length.out = 100)

disp.mod1.rss$df$D2.indiv.hr.175.P_x1_back <- (disp.scaled.range*sd(deer.data.disp$D2.indiv.hr.175.P)) + mean(deer.data.disp$D2.indiv.hr.175.P)

# plot prediction
disp.rss.plot <- ggplot(data = disp.mod1.rss$df, aes(x = exp(D2.indiv.hr.175.P_x1_back), y = log_rss)) +
       theme_bw() +
       geom_hline(yintercept = 0) +
       geom_vline(xintercept = exp(mean(deer.data.disp$D2.indiv.hr.175.P)), linetype = "dashed") +
       geom_line(size = 1.5, color = "black") +
       geom_ribbon(aes(ymin = lwr,
                       ymax = upr),
                   fill = "black",
                   alpha = 0.1) +
       theme(panel.grid = element_blank()) +
       ylab("") +
       xlab(expression(paste(D^2, " from indiv. HR locations (175-m)"))) +
       coord_cartesian(xlim = c(0, quantile(deer.data.disp$D2.indiv.hr.175, 0.99)),
                       ylim = c(-1.3, 1.6)) +
       scale_y_continuous(breaks = seq(-1, 1.5, 0.5)) +
       ggtitle("d")

#__________________________________________________________________________________________________
# 5. Multiple plots ----
#__________________________________________________________________________________________________

plot_grid(spring.rss.plot, summer.rss.plot, fall.rss.plot, disp.rss.plot,
          nrow = 1)

#__________________________________________________________________________________________________
# 6. Save image ----
#__________________________________________________________________________________________________

save.image(file = "SSF_predictions.RData")

#__________________________________________________________________________________________________
# 7. Individual variation ----
#__________________________________________________________________________________________________
# 7a. Spring ----
#__________________________________________________________________________________________________

# calculate responses for plotting
spring.iv.df <- data.frame()

for (x in 1:nrow(spring.coef$cond$burst_)) {
  
  burst.pred <- spring.coef$cond$burst_$D2.indiv.hr.175.P[x]*seq(min(deer.data.spring.s$D2.indiv.hr.175.P),
                                                                                              max(deer.data.spring.s$D2.indiv.hr.175.P),
                                                                                              length.out = 100)
  
  burst.df <- data.frame(range = seq(min(deer.data.spring.s$D2.indiv.hr.175.P),
                                     max(deer.data.spring.s$D2.indiv.hr.175.P),
                                     length.out = 100),
                         response = burst.pred,
                         group = rownames(spring.coef$cond$burst_)[x])
  
  spring.iv.df <- rbind(spring.iv.df, burst.df)
  
}

# back transform
spring.iv.df$range.back <- (spring.iv.df$range*sd(deer.data.spring$D2.indiv.hr.175.P)) + mean(deer.data.spring$D2.indiv.hr.175.P)

spring.iv.plot <- ggplot(data = spring.iv.df, aes(x = exp(range.back), y = response, group = group)) +
                     theme_bw() +
                     geom_hline(yintercept = 0) +
                     geom_line(color = "black",
                               alpha = 0.25) +
                     theme(panel.grid = element_blank(),
                           axis.title.y = element_blank()) +
                     geom_vline(xintercept = exp(mean(deer.data.spring$D2.indiv.hr.175.P)), linetype = "dashed") +
                     ylab("ln(Relative selection strength)") +
                     xlab(expression(paste(D^2, " from indiv. HR locations (175-m)"))) +
                     coord_cartesian(xlim = c(0, quantile(deer.data.spring$D2.indiv.hr.175, 0.99)),
                                     ylim = c(-1, 1)) +
                     ggtitle("a")

#__________________________________________________________________________________________________
# 7b. Summer ----
#__________________________________________________________________________________________________

# calculate responses for plotting
summer.iv.df <- data.frame()

for (x in 1:nrow(summer.coef$cond$burst_)) {
  
  burst.pred <- summer.coef$cond$burst_$D2.indiv.hr.175.P[x]*seq(min(deer.data.summer.s$D2.indiv.hr.175.P),
                                                                 max(deer.data.summer.s$D2.indiv.hr.175.P),
                                                                 length.out = 100)
  
  burst.df <- data.frame(range = seq(min(deer.data.summer.s$D2.indiv.hr.175.P),
                                     max(deer.data.summer.s$D2.indiv.hr.175.P),
                                     length.out = 100),
                         response = burst.pred,
                         group = rownames(summer.coef$cond$burst_)[x])
  
  summer.iv.df <- rbind(summer.iv.df, burst.df)
  
}

# back transform
summer.iv.df$range.back <- (summer.iv.df$range*sd(deer.data.summer$D2.indiv.hr.175.P)) + mean(deer.data.summer$D2.indiv.hr.175.P)

summer.iv.plot <- ggplot(data = summer.iv.df, aes(x = exp(range.back), y = response, group = group)) +
                     theme_bw() +
                     geom_hline(yintercept = 0) +
                     geom_vline(xintercept = exp(mean(deer.data.summer$D2.indiv.hr.175.P)), linetype = "dashed") +
                     geom_line(color = "black",
                               alpha = 0.25) +
                     theme(panel.grid = element_blank(),
                           axis.title.y = element_blank()) +
                     ylab("ln(Relative selection strength)") +
                     xlab(expression(paste(D^2, " from indiv. HR locations (175-m)"))) +
                     coord_cartesian(xlim = c(0, quantile(deer.data.summer$D2.indiv.hr.175, 0.99)),
                                     ylim = c(-0.6, 0.5)) +
                     ggtitle("b")

#__________________________________________________________________________________________________
# 7c. Fall ----
#__________________________________________________________________________________________________

# calculate responses for plotting
fall.iv.df <- data.frame()

for (x in 1:nrow(fall.coef$cond$burst_)) {
  
  burst.pred <- fall.coef$cond$burst_$D2.indiv.hr.350.P[x]*seq(min(deer.data.fall.s$D2.indiv.hr.350.P),
                                                                                              max(deer.data.fall.s$D2.indiv.hr.350.P),
                                                                                              length.out = 100)
  
  burst.df <- data.frame(range = seq(min(deer.data.fall.s$D2.indiv.hr.350.P),
                                     max(deer.data.fall.s$D2.indiv.hr.350.P),
                                     length.out = 100),
                         response = burst.pred,
                         group = rownames(fall.coef$cond$burst_)[x])
  
  fall.iv.df <- rbind(fall.iv.df, burst.df)
  
}

# back transform
fall.iv.df$range.back <- (fall.iv.df$range*sd(deer.data.fall$D2.indiv.hr.350.P)) + mean(deer.data.fall$D2.indiv.hr.350.P)

fall.iv.plot <- ggplot(data = fall.iv.df, aes(x = exp(range.back), y = response, group = group)) +
                     theme_bw() +
                     geom_hline(yintercept = 0) +
                     geom_vline(xintercept = exp(mean(deer.data.fall$D2.indiv.hr.350.P)), linetype = "dashed") +
                     geom_line(color = "black",
                               alpha = 0.25) +
                     theme(panel.grid = element_blank(),
                           axis.title.y = element_blank()) +
                     ylab("ln(Relative selection strength)") +
                     xlab(expression(paste(D^2, " from indiv. HR locations (350-m)"))) +
                     coord_cartesian(xlim = c(0, quantile(deer.data.fall$D2.indiv.hr.350, 0.99)),
                                     ylim = c(-1, 1)) +
                     scale_y_continuous(breaks = seq(-1, 1, 0.5)) +
                     ggtitle("c")

#__________________________________________________________________________________________________
# 7d. Disp ----
#__________________________________________________________________________________________________

# calculate responses for plotting
disp.iv.df <- data.frame()

for (x in 1:nrow(disp.coef$cond$burst_)) {
  
  burst.pred <- disp.coef$cond$burst_$D2.indiv.hr.175.P[x]*seq(min(deer.data.disp.s$D2.indiv.hr.175.P),
                                                                                              max(deer.data.disp.s$D2.indiv.hr.175.P),
                                                                                              length.out = 100)
  
  burst.df <- data.frame(range = seq(min(deer.data.disp.s$D2.indiv.hr.175.P),
                                     max(deer.data.disp.s$D2.indiv.hr.175.P),
                                     length.out = 100),
                         response = burst.pred,
                         group = rownames(disp.coef$cond$burst_)[x])
  
  disp.iv.df <- rbind(disp.iv.df, burst.df)
  
}

# back transform
disp.iv.df$range.back <- (disp.iv.df$range*sd(deer.data.disp$D2.indiv.hr.175.P)) + mean(deer.data.disp$D2.indiv.hr.175.P)

disp.iv.plot <- ggplot(data = disp.iv.df, aes(x = exp(range.back), y = response, group = group)) +
                     theme_bw() +
                     geom_hline(yintercept = 0) +
                     geom_vline(xintercept = exp(mean(deer.data.disp$D2.indiv.hr.175.P)), linetype = "dashed") +
                     geom_line(color = "black",
                               alpha = 0.25) +
                     theme(panel.grid = element_blank(),
                           axis.title.y = element_blank()) +
                     ylab("ln(Relative selection strength)") +
                     xlab(expression(paste(D^2, " from indiv. HR locations (175-m)"))) +
                     coord_cartesian(xlim = c(0, quantile(deer.data.disp$D2.indiv.hr.175, 0.99)),
                                     ylim = c(-1.3, 1.6)) +
                     scale_y_continuous(breaks = seq(-1, 1.5, 0.5)) +
                     ggtitle("d")

cowplot::plot_grid(spring.iv.plot, summer.iv.plot, fall.iv.plot, disp.iv.plot,
                   nrow = 1)

#__________________________________________________________________________________________________
# 8. Save image ----
#__________________________________________________________________________________________________

save.image(file = "SSF_predictions.RData")
