# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 10 - PCA plot
# Author: Nathan D. Hooven
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Email: nathan.hooven@uky.edu
# Date started: 10 May 2022
# Date completed: 27 May 2022
# Date modified: 
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(factoextra)
library(cowplot)

#__________________________________________________________________________________________________
# 2. Read in data ----
#__________________________________________________________________________________________________

# load
load("SSF_models.RData")

# create "burst_" column
all.pre.data$burst_ <- paste0(all.pre.data$Deer, ".", all.pre.data$MV.ID)

# merge Season column from deer.data.2
deer.data.3 <- deer.data.2 %>% dplyr::select(burst_, Season)

# for loop because left_join runs out of memory
all.pre.data.1 <- data.frame()

for (x in unique(all.pre.data$burst_)) {
  
  indiv.burst.pre <- all.pre.data %>% filter(burst_ == x)
  
  indiv.burst.mv <- deer.data.3[deer.data.3$burst_ == x,][1, 2]
  
  indiv.burst.pre$Season <- indiv.burst.mv
  
  # bind to master df
  all.pre.data.1 <- rbind(all.pre.data.1, indiv.burst.pre)
  
}

#__________________________________________________________________________________________________
# 3. Fit PCAs ----
#__________________________________________________________________________________________________
# 3a. 175-m ---
#__________________________________________________________________________________________________

pre.PCA.175 <- prcomp(all.pre.data.1[ ,c(17:19, 21:24)], scale = TRUE)

# variance explained by PCs
get_eig(pre.PCA.175)

# bind first two PCs to df
all.pre.data.1$PC1.175 <- get_pca_ind(pre.PCA.175)$coord[ ,1]
all.pre.data.1$PC2.175 <- get_pca_ind(pre.PCA.175)$coord[ ,2]

#__________________________________________________________________________________________________
# 3b. 250-m ---
#__________________________________________________________________________________________________

pre.PCA.250 <- prcomp(all.pre.data.1[ ,c(25:27, 29:32)], scale = TRUE)

# variance explained by PCs
get_eig(pre.PCA.250)

# bind first two PCs to df
all.pre.data.1$PC1.250 <- get_pca_ind(pre.PCA.250)$coord[ ,1]
all.pre.data.1$PC2.250 <- get_pca_ind(pre.PCA.250)$coord[ ,2]

#__________________________________________________________________________________________________
# 3c. 350-m ---
#__________________________________________________________________________________________________

pre.PCA.350 <- prcomp(all.pre.data.1[ ,c(33:35, 37:40)], scale = TRUE)

# variance explained by PCs
get_eig(pre.PCA.350)

# bind first two PCs to df
all.pre.data.1$PC1.350 <- get_pca_ind(pre.PCA.350)$coord[ ,1]
all.pre.data.1$PC2.350 <- get_pca_ind(pre.PCA.350)$coord[ ,2]

#__________________________________________________________________________________________________
# 3d. 500-m ---
#__________________________________________________________________________________________________

pre.PCA.500 <- prcomp(all.pre.data.1[ ,c(41:43, 45:48)], scale = TRUE)

# variance explained by PCs
get_eig(pre.PCA.500)

# bind first two PCs to df
all.pre.data.1$PC1.500 <- get_pca_ind(pre.PCA.500)$coord[ ,1]
all.pre.data.1$PC2.500 <- get_pca_ind(pre.PCA.500)$coord[ ,2]

#__________________________________________________________________________________________________
# 4. Individual specialization by centroid  ----
#__________________________________________________________________________________________________

all.pre.data.summary <- all.pre.data.1 %>% group_by(burst_) %>%
                                           summarize_at(c(17:19, 21:24,
                                                       25:27, 29:32,
                                                       33:35, 37:40,
                                                       41:43, 45:48), 
                                                     mean) %>% 
                                           left_join(distinct(deer.data.3), by = "burst_")

# convert Season to factor
all.pre.data.summary$Season <- factor(all.pre.data.summary$Season,
                                      levels = c("Spring", "Summer", "Fall", "Disp"),
                                      labels = c("Spring", "Summer", "Fall", "Dispersal"))

# predict using PCAs
all.pre.data.summary$PC1.175 <- predict(pre.PCA.175, all.pre.data.summary)[ ,1]
all.pre.data.summary$PC2.175 <- predict(pre.PCA.175, all.pre.data.summary)[ ,2]

all.pre.data.summary$PC1.250 <- predict(pre.PCA.250, all.pre.data.summary)[ ,1]
all.pre.data.summary$PC2.250 <- predict(pre.PCA.250, all.pre.data.summary)[ ,2]

all.pre.data.summary$PC1.350 <- predict(pre.PCA.350, all.pre.data.summary)[ ,1]
all.pre.data.summary$PC2.350 <- predict(pre.PCA.350, all.pre.data.summary)[ ,2]

all.pre.data.summary$PC1.500 <- predict(pre.PCA.500, all.pre.data.summary)[ ,1]
all.pre.data.summary$PC2.500 <- predict(pre.PCA.500, all.pre.data.summary)[ ,2]

# overall HR means
all.pre.data.summary.overall <- as.data.frame(t(as.matrix(colMeans(all.pre.data.summary[ ,2:29]))))

all.pre.data.summary.overall$PC1.175 <- predict(pre.PCA.175, all.pre.data.summary.overall)[ ,1]
all.pre.data.summary.overall$PC2.175 <- predict(pre.PCA.175, all.pre.data.summary.overall)[ ,2]

all.pre.data.summary.overall$PC1.250 <- predict(pre.PCA.250, all.pre.data.summary.overall)[ ,1]
all.pre.data.summary.overall$PC2.250 <- predict(pre.PCA.250, all.pre.data.summary.overall)[ ,2]

all.pre.data.summary.overall$PC1.350 <- predict(pre.PCA.350, all.pre.data.summary.overall)[ ,1]
all.pre.data.summary.overall$PC2.350 <- predict(pre.PCA.350, all.pre.data.summary.overall)[ ,2]

all.pre.data.summary.overall$PC1.500 <- predict(pre.PCA.500, all.pre.data.summary.overall)[ ,1]
all.pre.data.summary.overall$PC2.500 <- predict(pre.PCA.500, all.pre.data.summary.overall)[ ,2]


# plot altogether
# 175-m
pca.plot.175 <- ggplot(data = all.pre.data.summary, aes(x = PC1.175, y = PC2.175)) +
                theme_bw() +
                geom_vline(xintercept = 0) +
                geom_hline(yintercept = 0) +
                geom_point(aes(fill = Season,
                               shape = Season),
                           size = 3) +
                scale_fill_manual(values = c("darkgreen", "orange", "red", "black")) +
                scale_shape_manual(values = c(21, 22, 24, 23)) +
                geom_point(data = all.pre.data.summary.overall,
                           aes(x = PC1.175, y = PC2.175),
                           size = 4,
                           shape = 21,
                           fill = "white",
                           stroke = 2) +
                theme(legend.position = "none", 
                      panel.grid = element_blank()) +
                xlab("PC1 (42.9%)") +
                ylab("PC2 (16.5%)") +
                ggtitle("a")

# 250-m
pca.plot.250 <- ggplot(data = all.pre.data.summary, aes(x = PC1.250, y = PC2.250)) +
                theme_bw() +
                geom_vline(xintercept = 0) +
                geom_hline(yintercept = 0) +
                geom_point(aes(fill = Season,
                               shape = Season),
                           size = 3) +
                scale_fill_manual(values = c("darkgreen", "orange", "red", "black")) +
                scale_shape_manual(values = c(21, 22, 24, 23)) +
                geom_point(data = all.pre.data.summary.overall,
                           aes(x = PC1.250, y = PC2.250),
                           size = 4,
                           shape = 21,
                           fill = "white",
                           stroke = 2) +
                theme(legend.position = "none", 
                      panel.grid = element_blank()) +
                xlab("PC1 (38.7%)") +
                ylab("PC2 (18.7%)") +
                ggtitle("b")

# 350-m
pca.plot.350 <- ggplot(data = all.pre.data.summary, aes(x = PC1.350, y = PC2.350)) +
                theme_bw() +
                geom_vline(xintercept = 0) +
                geom_hline(yintercept = 0) +
                geom_point(aes(fill = Season,
                               shape = Season),
                           size = 3) +
                scale_fill_manual(values = c("darkgreen", "orange", "red", "black")) +
                scale_shape_manual(values = c(21, 22, 24, 23)) +
                geom_point(data = all.pre.data.summary.overall,
                           aes(x = PC1.350, y = PC2.350),
                           size = 4,
                           shape = 21,
                           fill = "white",
                           stroke = 2) +
                theme(legend.position = "none", 
                      panel.grid = element_blank()) +
                xlab("PC1 (38.6%)") +
                ylab("PC2 (18.7%)") +
                ggtitle("c")

# 500-m
pca.plot.500 <- ggplot(data = all.pre.data.summary, aes(x = PC1.500, y = PC2.500)) +
                theme_bw() +
                geom_vline(xintercept = 0) +
                geom_hline(yintercept = 0) +
                geom_point(aes(fill = Season,
                               shape = Season),
                           size = 3) +
                scale_fill_manual(values = c("darkgreen", "orange", "red", "black")) +
                scale_shape_manual(values = c(21, 22, 24, 23)) +
                geom_point(data = all.pre.data.summary.overall,
                           aes(x = PC1.500, y = PC2.500),
                           size = 4,
                           shape = 21,
                           fill = "white",
                           stroke = 2) +
                theme(legend.position = "none", 
                      panel.grid = element_blank()) +
                xlab("PC1 (42.3%)") +
                ylab("PC2 (17.2%)") +
                ggtitle("d") +
                scale_y_continuous(breaks = c(-1, 0, 1, 2),
                                   labels = c("-1.0", "0.0", "1.0", "2.0"))

# plot altogether
plot_grid(pca.plot.175, pca.plot.250, pca.plot.350, pca.plot.500, nrow = 2)

save.image("PCA_plots.RData")

#__________________________________________________________________________________________________
# 5. Convex hull visualization ----

# correct factor order
all.pre.data.1$Season <- factor(all.pre.data.1$Season,
                                levels = c("Spring", "Summer", "Fall", "Disp"))

#__________________________________________________________________________________________________
# 5a. 175-m ----
#__________________________________________________________________________________________________

# create dataset
hull.175 <- all.pre.data.1 %>% group_by(burst_) %>%
                               slice(chull(PC1.175, PC2.175))

# spring
hull.plot.175.spring <- ggplot(data = hull.175, aes(x = PC1.175, y = PC2.175)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.175[hull.175$Season != "Spring", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.175[hull.175$Season == "Spring", ],
                                            color = "darkgreen",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (42.9%)") +
                               ylab("PC2 (16.5%)") +
                               scale_y_continuous(breaks = c(-2, 0, 2, 4))

# summer
hull.plot.175.summer <- ggplot(data = hull.175, aes(x = PC1.175, y = PC2.175)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.175[hull.175$Season != "Summer", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.175[hull.175$Season == "Summer", ],
                                            color = "orange",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (42.9%)") +
                               ylab("PC2 (16.5%)") +
                               scale_y_continuous(breaks = c(-2, 0, 2, 4))

# fall
hull.plot.175.fall <- ggplot(data = hull.175, aes(x = PC1.175, y = PC2.175)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.175[hull.175$Season != "Fall", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.175[hull.175$Season == "Fall", ],
                                            color = "red",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (42.9%)") +
                               ylab("PC2 (16.5%)") +
                               scale_y_continuous(breaks = c(-2, 0, 2, 4))

# disp
hull.plot.175.disp <- ggplot(data = hull.175, aes(x = PC1.175, y = PC2.175)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.175[hull.175$Season != "Disp", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.175[hull.175$Season == "Disp", ],
                                            color = "black",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (42.9%)") +
                               ylab("PC2 (16.5%)") +
                               scale_y_continuous(breaks = c(-2, 0, 2, 4))

#__________________________________________________________________________________________________
# 5b. 250-m ----
#__________________________________________________________________________________________________

# create dataset
hull.250 <- all.pre.data.1 %>% group_by(burst_) %>%
                               slice(chull(PC1.250, PC2.250))

# spring
hull.plot.250.spring <- ggplot(data = hull.250, aes(x = PC1.250, y = PC2.250)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.250[hull.250$Season != "Spring", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.250[hull.250$Season == "Spring", ],
                                            color = "darkgreen",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (38.7%)") +
                               ylab("PC2 (18.7%)")

# summer
hull.plot.250.summer <- ggplot(data = hull.250, aes(x = PC1.250, y = PC2.250)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.250[hull.250$Season != "Summer", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.250[hull.250$Season == "Summer", ],
                                            color = "orange",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (38.7%)") +
                               ylab("PC2 (18.7%)")

# fall
hull.plot.250.fall <- ggplot(data = hull.250, aes(x = PC1.250, y = PC2.250)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.250[hull.250$Season != "Fall", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.250[hull.250$Season == "Fall", ],
                                            color = "red",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (38.7%)") +
                               ylab("PC2 (18.7%)")

# disp
hull.plot.250.disp <- ggplot(data = hull.250, aes(x = PC1.250, y = PC2.250)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.250[hull.250$Season != "Disp", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.250[hull.250$Season == "Disp", ],
                                            color = "black",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (38.7%)") +
                               ylab("PC2 (18.7%)")

#__________________________________________________________________________________________________
# 5c. 350-m ----
#__________________________________________________________________________________________________

# create dataset
hull.350 <- all.pre.data.1 %>% group_by(burst_) %>%
                               slice(chull(PC1.350, PC2.350))

# spring
hull.plot.350.spring <- ggplot(data = hull.350, aes(x = PC1.350, y = PC2.350)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.350[hull.350$Season != "Spring", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.350[hull.350$Season == "Spring", ],
                                            color = "darkgreen",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (38.6%)") +
                               ylab("PC2 (18.7%)")

# summer
hull.plot.350.summer <- ggplot(data = hull.350, aes(x = PC1.350, y = PC2.350)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.350[hull.350$Season != "Summer", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.350[hull.350$Season == "Summer", ],
                                            color = "orange",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (38.6%)") +
                               ylab("PC2 (18.7%)")

# fall
hull.plot.350.fall <- ggplot(data = hull.350, aes(x = PC1.350, y = PC2.350)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.350[hull.350$Season != "Fall", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.350[hull.350$Season == "Fall", ],
                                            color = "red",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (38.6%)") +
                               ylab("PC2 (18.7%)")

# disp
hull.plot.350.disp <- ggplot(data = hull.350, aes(x = PC1.350, y = PC2.350)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.350[hull.350$Season != "Disp", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.350[hull.350$Season == "Disp", ],
                                            color = "black",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (38.6%)") +
                               ylab("PC2 (18.7%)")

#__________________________________________________________________________________________________
# 5d. 500-m ----
#__________________________________________________________________________________________________

# create dataset
hull.500 <- all.pre.data.1 %>% group_by(burst_) %>%
                               slice(chull(PC1.500, PC2.500))

# spring
hull.plot.500.spring <- ggplot(data = hull.500, aes(x = PC1.500, y = PC2.500)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.500[hull.500$Season != "Spring", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.500[hull.500$Season == "Spring", ],
                                            color = "darkgreen",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (42.3%)") +
                               ylab("PC2 (17.2%)")

# summer
hull.plot.500.summer <- ggplot(data = hull.500, aes(x = PC1.500, y = PC2.500)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.500[hull.500$Season != "Summer", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.500[hull.500$Season == "Summer", ],
                                            color = "orange",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (42.3%)") +
                               ylab("PC2 (17.2%)")

# fall
hull.plot.500.fall <- ggplot(data = hull.500, aes(x = PC1.500, y = PC2.500)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.500[hull.500$Season != "Fall", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.500[hull.500$Season == "Fall", ],
                                            color = "red",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (42.3%)") +
                               ylab("PC2 (17.2%)")

# disp
hull.plot.500.disp <- ggplot(data = hull.500, aes(x = PC1.500, y = PC2.500)) +
                               theme_bw() +
                               geom_vline(xintercept = 0) +
                               geom_hline(yintercept = 0) +
                               geom_polygon(data = hull.500[hull.500$Season != "Disp", ],
                                            color = "lightgray",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               geom_polygon(data = hull.500[hull.500$Season == "Disp", ],
                                            color = "black",
                                            size = 0.5,
                                            fill = NA,
                                            aes(group = burst_)) +
                               theme(legend.position = "none", 
                                     panel.grid = element_blank()) +
                               xlab("PC1 (42.3%)") +
                               ylab("PC2 (17.2%)")

#__________________________________________________________________________________________________
# 5e. Plot altogether ----
#__________________________________________________________________________________________________

plot_grid(hull.plot.175.spring, hull.plot.175.summer, hull.plot.175.fall, hull.plot.175.disp,
          hull.plot.250.spring, hull.plot.250.summer, hull.plot.250.fall, hull.plot.250.disp,
          hull.plot.350.spring, hull.plot.350.summer, hull.plot.350.fall, hull.plot.350.disp,
          hull.plot.500.spring, hull.plot.500.summer, hull.plot.500.fall, hull.plot.500.disp,
          nrow = 4)
