# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 10 - Prediction plots
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
library(lubridate)     # work with dates
library(AICcmodavg)    # model selection
library(amt)           # calculate RSS
library(cowplot)       # multiple plots
library(glmmTMB)       # mixed models

#__________________________________________________________________________________________________
# 2. Read in data ----
#__________________________________________________________________________________________________

load("SSF_models_NEW.RData")

#__________________________________________________________________________________________________
# 3. Relative selection strength prediction plots ----
#__________________________________________________________________________________________________
# 3a. Spring ----
#__________________________________________________________________________________________________

# variance-covariance matrix
vcov.spring <- vcov(spring.models[[42]])

m.vcov.spring <- as.matrix(vcov.spring$cond)

m.vcov.spring <- m.vcov.spring[-1, -1]         # remove intercept

# create a sequence for the variable of interest
spring.seq <- seq(min(deer.data.spring.s$D2.indiv.hr.175.P), max(deer.data.spring.s$D2.indiv.hr.175.P), length.out = 100)

spring.x2 <- data.frame(dForest.P = 0,
                        dStream = 0,
                        dStream.Q = 0,
                        sl_ = 0,
                        D2.all.hr.175.P = 0,
                        D2.indiv.hr.175.P = 0)

spring.x1 <- data.frame(dForest.P = 0,
                        dStream = 0,
                        dStream.Q = 0,
                        sl_ = 0,
                        D2.all.hr.175.P = 0,
                        D2.indiv.hr.175.P = spring.seq)

# model matrices
# x2 (avg)
mm.spring.x2 <- model.matrix(~ dForest.P +
                               dStream +
                               dStream.Q +
                               sl_ +
                               D2.all.hr.175.P +
                               D2.indiv.hr.175.P,
                             spring.x2)


# x1 (varying)
mm.spring.x1 <- model.matrix(~ dForest.P +
                                 dStream +
                                 dStream.Q +
                                 sl_ +
                                 D2.all.hr.175.P +
                                 D2.indiv.hr.175.P,
                               spring.x1)

# calculate RSS for base + all + indiv model
spring.x1$rss <-  spring.seq*spring.models[[42]]$fit$par[6]

# create confidence intervals
delta.mm.spring.X1 <- sweep(data.matrix(mm.spring.x1), 2, data.matrix(mm.spring.x2))
  
spring.vars.X1 <- delta.mm.spring.X1[ ,1:ncol(spring.x2) + 1] %*% m.vcov.spring %*% t(delta.mm.spring.X1[ ,1:ncol(spring.x2) + 1])
  
# create standard errors
spring.sds.X1 <- sqrt(diag(spring.vars.X1))
  
# create 95% cIs
spring.x1$low <- spring.x1$rss - 1.96 * spring.sds.X1
spring.x1$upp <- spring.x1$rss + 1.96 * spring.sds.X1

# model average
spring.x1[ ,7:9] <- spring.x1[ ,7:9] * 0.05

# determine range for back-transformation
spring.scaled.range <- seq(min(deer.data.spring.s$D2.indiv.hr.175.P),
                           max(deer.data.spring.s$D2.indiv.hr.175.P),
                           length.out = 100)

spring.x1$x.back <- (spring.scaled.range*sd(deer.data.spring$D2.indiv.hr.175.P)) + mean(deer.data.spring$D2.indiv.hr.175.P)

# plot prediction
spring.rss.plot <- ggplot(data = spring.x1, aes(x = exp(x.back), 
                                                y = rss)) +
       theme_bw() +
       geom_hline(yintercept = 0) +
       geom_vline(xintercept = exp(mean(deer.data.spring$D2.indiv.hr.175.P)), linetype = "dashed") +
       geom_line(size = 1.5, color = "darkgreen") +
       geom_ribbon(aes(ymin = low,
                       ymax = upp),
                   fill = "darkgreen",
                   alpha = 0.1) +
       theme(panel.grid = element_blank()) +
       ylab("ln(Relative selection strength)") +
       xlab(expression(paste(D^2, " from indiv. HR locations (175-m)"))) +
       coord_cartesian(xlim = c(0, quantile(deer.data.spring$D2.indiv.hr.175, 0.99)),
                       ylim = c(-0.01, 0.015)) +
       ggtitle("a")

#__________________________________________________________________________________________________
# 3b. Summer ----
#__________________________________________________________________________________________________

# variance-covariance matrix
vcov.summer <- vcov(summer.models[[42]])

m.vcov.summer <- as.matrix(vcov.summer$cond)

m.vcov.summer <- m.vcov.summer[-1, -1]         # remove intercept

# create a sequence for the variable of interest
summer.seq <- seq(min(deer.data.summer.s$D2.indiv.hr.175.P), max(deer.data.summer.s$D2.indiv.hr.175.P), length.out = 100)

summer.x2 <- data.frame(dForest.P = 0,
                        dStream = 0,
                        sl_ = 0,
                        D2.all.hr.250.P = 0,
                        D2.indiv.hr.175.P = 0)

summer.x1 <- data.frame(dForest.P = 0,
                        dStream = 0,
                        sl_ = 0,
                        D2.all.hr.250.P = 0,
                        D2.indiv.hr.175.P = summer.seq)

# model matrices
# x2 (avg)
mm.summer.x2 <- model.matrix(~ dForest.P +
                               dStream +
                               sl_ +
                               D2.all.hr.250.P +
                               D2.indiv.hr.175.P,
                             summer.x2)


# x1 (varying)
mm.summer.x1 <- model.matrix(~ dForest.P +
                                 dStream +
                                 sl_ +
                                 D2.all.hr.250.P +
                                 D2.indiv.hr.175.P,
                               summer.x1)

# calculate RSS for base + all + indiv model
summer.x1$rss <-  summer.seq*summer.models[[42]]$fit$par[5]

# create confidence intervals
delta.mm.summer.X1 <- sweep(data.matrix(mm.summer.x1), 2, data.matrix(mm.summer.x2))
  
summer.vars.X1 <- delta.mm.summer.X1[ ,1:ncol(summer.x2) + 1] %*% m.vcov.summer %*% t(delta.mm.summer.X1[ ,1:ncol(summer.x2) + 1])
  
# create standard errors
summer.sds.X1 <- sqrt(diag(summer.vars.X1))
  
# create 95% cIs
summer.x1$low <- summer.x1$rss - 1.96 * summer.sds.X1
summer.x1$upp <- summer.x1$rss + 1.96 * summer.sds.X1

# model average
summer.x1[ ,6:8] <- summer.x1[ ,6:8] * 0.23

# determine range for back-transformation
summer.scaled.range <- seq(min(deer.data.summer.s$D2.indiv.hr.175.P),
                           max(deer.data.summer.s$D2.indiv.hr.175.P),
                           length.out = 100)

summer.x1$x.back <- (summer.scaled.range*sd(deer.data.summer$D2.indiv.hr.175.P)) + mean(deer.data.summer$D2.indiv.hr.175.P)

# plot prediction
summer.rss.plot <- ggplot(data = summer.x1, aes(x = exp(x.back), y = rss)) +
       theme_bw() +
       geom_hline(yintercept = 0) +
       geom_vline(xintercept = exp(mean(deer.data.summer$D2.indiv.hr.175.P)), linetype = "dashed") +
       geom_line(size = 1.5, color = "orange") +
       geom_ribbon(aes(ymin = low,
                       ymax = upp),
                   fill = "orange",
                   alpha = 0.1) +
       theme(panel.grid = element_blank()) +
       ylab("") +
       xlab(expression(paste(D^2, " from indiv. HR locations (175-m)"))) +
       coord_cartesian(xlim = c(0, quantile(deer.data.summer$D2.indiv.hr.175, 0.99)),
                       ylim = c(-0.21, 0.15)) +
       ggtitle("b")

#__________________________________________________________________________________________________
# 3c. Fall ----
#__________________________________________________________________________________________________

# variance-covariance matrix
vcov.fall <- vcov(fall.models[[42]])

m.vcov.fall <- as.matrix(vcov.fall$cond)

m.vcov.fall <- m.vcov.fall[-1, -1]         # remove intercept

# create a sequence for the variable of interest
fall.seq <- seq(min(deer.data.fall.s$D2.indiv.hr.350.P), max(deer.data.fall.s$D2.indiv.hr.350.P), length.out = 100)

fall.x2 <- data.frame(dForest.P = 0,
                        dStream = 0,
                        dAg.P = 0,
                        dRoad.P = 0,
                        sl_ = 0,
                        D2.all.hr.175.P = 0,
                        D2.indiv.hr.350.P = 0)

fall.x1 <- data.frame(dForest.P = 0,
                        dStream = 0,
                        dAg.P = 0,
                        dRoad.P = 0,
                        sl_ = 0,
                        D2.all.hr.175.P = 0,
                        D2.indiv.hr.350.P = fall.seq)

# model matrices
# x2 (avg)
mm.fall.x2 <- model.matrix(~ dForest.P +
                               dStream +
                             dAg.P +
                             dRoad.P +
                               sl_ +
                               D2.all.hr.175.P +
                               D2.indiv.hr.350.P,
                             fall.x2)

# x1 (varying)
mm.fall.x1 <- model.matrix(~ dForest.P +
                               dStream +
                             dAg.P +
                             dRoad.P +
                               sl_ +
                               D2.all.hr.175.P +
                               D2.indiv.hr.350.P,
                               fall.x1)

# calculate RSS for base + all + indiv model
fall.x1$rss <-  fall.seq*fall.models[[42]]$fit$par[7]

# create confidence intervals
delta.mm.fall.X1 <- sweep(data.matrix(mm.fall.x1), 2, data.matrix(mm.fall.x2))
  
fall.vars.X1 <- delta.mm.fall.X1[ ,1:ncol(fall.x2) + 1] %*% m.vcov.fall %*% t(delta.mm.fall.X1[ ,1:ncol(fall.x2) + 1])
  
# create standard errors
fall.sds.X1 <- sqrt(diag(fall.vars.X1))
  
# create 95% cIs
fall.x1$low <- fall.x1$rss - 1.96 * fall.sds.X1
fall.x1$upp <- fall.x1$rss + 1.96 * fall.sds.X1

# model average
fall.x1[ ,8:10] <- fall.x1[ ,8:10] * 0.31

# determine range for back-transformation
fall.scaled.range <- seq(min(deer.data.fall.s$D2.indiv.hr.350.P),
                           max(deer.data.fall.s$D2.indiv.hr.350.P),
                           length.out = 100)

fall.x1$x.back <- (fall.scaled.range*sd(deer.data.fall$D2.indiv.hr.350.P)) + mean(deer.data.fall$D2.indiv.hr.350.P)

# plot prediction
fall.rss.plot <- ggplot(data = fall.x1, aes(x = exp(x.back), y = rss)) +
       theme_bw() +
       geom_hline(yintercept = 0) +
       geom_vline(xintercept = 35, linetype = "dashed") +
       geom_line(size = 1.5, color = "red") +
       geom_ribbon(aes(ymin = low,
                       ymax = upp),
                   fill = "red",
                   alpha = 0.1) +
       theme(panel.grid = element_blank()) +
       ylab("") +
       xlab(expression(paste(D^2, " from indiv. HR locations (350-m)"))) +
       coord_cartesian(xlim = c(0, quantile(deer.data.fall$D2.indiv.hr.350, 0.99)),
                       ylim = c(-0.15, 0.10)) +
       scale_y_continuous(breaks = seq(-0.20, 0.10, 0.05)) +
       ggtitle("c")

#__________________________________________________________________________________________________
# 3d. Dispersal ----
#__________________________________________________________________________________________________

# variance-covariance matrix
vcov.disp <- vcov(disp.models[[42]])

m.vcov.disp <- as.matrix(vcov.disp$cond)

m.vcov.disp <- m.vcov.disp[-1, -1]         # remove intercept

# create a sequence for the variable of interest
disp.seq <- seq(min(deer.data.disp.s$D2.indiv.hr.175.P), max(deer.data.disp.s$D2.indiv.hr.175.P), length.out = 100)

disp.x2 <- data.frame(dForest.P = 0,
                        dStream.P = 0,
                        sl_ = 0,
                        D2.all.hr.175 = 0,
                        D2.indiv.hr.175.P = 0)

disp.x1 <- data.frame(dForest.P = 0,
                        dStream.P = 0,
                        sl_ = 0,
                        D2.all.hr.175 = 0,
                        D2.indiv.hr.175.P = disp.seq)

# model matrices
# x2 (avg)
mm.disp.x2 <- model.matrix(~ dForest.P +
                               dStream.P +
                               sl_ +
                               D2.all.hr.175 +
                               D2.indiv.hr.175.P,
                             disp.x2)

# x1 (varying)
mm.disp.x1 <- model.matrix(~ dForest.P +
                               dStream.P +
                               sl_ +
                               D2.all.hr.175 +
                               D2.indiv.hr.175.P,
                               disp.x1)

# calculate RSS for base + all + indiv model
disp.x1$rss <-  disp.seq*disp.models[[42]]$fit$par[5]

# create confidence intervals
delta.mm.disp.X1 <- sweep(data.matrix(mm.disp.x1), 2, data.matrix(mm.disp.x2))
  
disp.vars.X1 <- delta.mm.disp.X1[ ,1:ncol(disp.x2) + 1] %*% m.vcov.disp %*% t(delta.mm.disp.X1[ ,1:ncol(disp.x2) + 1])
  
# create standard errors
disp.sds.X1 <- sqrt(diag(disp.vars.X1))
  
# create 95% cIs
disp.x1$low <- disp.x1$rss - 1.96 * disp.sds.X1
disp.x1$upp <- disp.x1$rss + 1.96 * disp.sds.X1

# model average
disp.x1[ ,6:8] <- disp.x1[ ,6:8] * 1.00

# determine range for back-transformation
disp.scaled.range <- seq(min(deer.data.disp.s$D2.indiv.hr.175.P),
                           max(deer.data.disp.s$D2.indiv.hr.175.P),
                           length.out = 100)

disp.x1$x.back <- (disp.scaled.range*sd(deer.data.disp$D2.indiv.hr.175.P)) + mean(deer.data.disp$D2.indiv.hr.175.P)

# plot prediction
disp.rss.plot <- ggplot(data = disp.x1, aes(x = exp(x.back), y = rss)) +
       theme_bw() +
       geom_hline(yintercept = 0) +
       geom_vline(xintercept = exp(mean(deer.data.disp$D2.indiv.hr.175.P)), linetype = "dashed") +
       geom_line(size = 1.5, color = "black") +
       geom_ribbon(aes(ymin = low,
                       ymax = upp),
                   fill = "black",
                   alpha = 0.1) +
       theme(panel.grid = element_blank()) +
       ylab("") +
       xlab(expression(paste(D^2, " from indiv. HR locations (175-m)"))) +
       coord_cartesian(xlim = c(0, quantile(deer.data.disp$D2.indiv.hr.175, 0.99)),
                       ylim = c(-1, 1)) +
       scale_y_continuous(breaks = seq(-1.0, 1.0, 0.5)) +
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
