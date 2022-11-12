# Title: Natal habitat preference during EHRM in WTD
# Subtitle: 1 - Raw data processing
# Author: Nathan D. Hooven
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Email: nathan.hooven@uky.edu
# Date started: 5 Dec 2020
# Date completed: 5 Dec 2020
# Date modified: 
# R version: 3.6.2

#__________________________________________________________________________________________________
# 1. Load required packages ----
#__________________________________________________________________________________________________

library(tidyverse)     # manage data
library(sp)            # work with spatial data
library(amt)           # create random steps
library(adehabitatHR)  # create home ranges
library(adehabitatLT)  # create ltraj objects
library(rgdal)         # save shapefiles

#__________________________________________________________________________________________________
# 2. Read in raw data ----
#__________________________________________________________________________________________________

rawdata <- read.csv("alldeerlocs.csv")

# select only columns we need
deerdata <- rawdata %>% select(Deer, Sex, Age, Date, Hour, Minute, Latitude, Longitude, HDOP, X2D_3D, UTME, UTMN)

#__________________________________________________________________________________________________
# 3. Process data ----
#__________________________________________________________________________________________________

# remove blank rows (by dropping all rows with NAs for UTME, an essential variable)
deerdata <- deerdata %>% drop_na(UTME)

# keep only relocations with HDOP <5 if 2-dimensional fix and HDOP <6 if 3-dimensional fix
deerdata <- deerdata %>% filter((X2D_3D == 2 & HDOP < 5) | (X2D_3D == 3 & HDOP < 6))

# create timestamp variable
deerdata$Minute <- paste0("0", deerdata$Minute)

deerdata <- deerdata %>% mutate(Timestamp = paste0(Date, " ", Hour, ":", Minute))

deerdata$Timestamp <- as.POSIXct(deerdata$Timestamp, format = "%m/%d/%Y %H:%M")

# Drop NAs
deerdata <- deerdata %>% drop_na(Deer)

# replace 247_1 with 247 and 252_1 with 252
deerdata$Deer[which(deerdata$Deer == "247_1")] <- "247"
deerdata$Deer[which(deerdata$Deer == "252_1")] <- "252"

# create final dataset
final.data <- deerdata %>% select(Deer, Sex, Age, Longitude, Latitude, UTME, UTMN, Timestamp)

# write to .csv
write.csv(final.data, "deerdata_1.csv")

#__________________________________________________________________________________________________
# 4. Remove post-capture relocations (7 days) ----
#__________________________________________________________________________________________________

# create list of deer names
deer.names <- as.vector(unique(final.data$Deer))

# create new data frame
final.data1 <- data.frame()

# for loop to remove post-capture relocations for each deer
for (i in deer.names){
  
  DeerID <- i
  
  only.data <- final.data %>% filter(Deer == DeerID)
  
  # determine first relocation
  first.relocation <- only.data[1, 8]
  
  # add 7 days
  end.relocation <- first.relocation + 7*24*60*60
  
  # only keep relocations after end.relocation
  only.data <- only.data %>% filter(Timestamp > end.relocation)
  
  # bind to final.data1
  final.data1 <- rbind(final.data1, only.data)
  
}

#__________________________________________________________________________________________________
# 5. Determine other erroneous relocations ----
#__________________________________________________________________________________________________

# for loop to fit an ltraj to each deer's relocations
for (i in deer.names){
  
  DeerID <- i
  
  only.data <- final.data1 %>% filter(Deer == DeerID)
  
  only.data <- final.data1[final.data1$Deer == Deer,]
  
  only.data <- only.data %>% drop_na(Deer)
  
  # only keep unique rows to prevent ltraj from throwing an error
  only.data <- only.data %>% distinct(Timestamp, .keep_all = TRUE)
  
  deer.ltraj <- as.ltraj(xy = only.data[c(6,7)], 
                         date = only.data$Timestamp,
                         id = Deer,
                         burst = nrow(only.data),
                         proj4string = CRS("+proj=utm +zone=16 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
  
  plot(deer.ltraj, main = Deer)
  
}

# individual relocation censoring

final.data2 <- data.frame()

# 157 - no censoring
only.data <- final.data1 %>% filter(Deer == "157")
only.data <- only.data %>% drop_na(Deer)
final.data2 <- rbind(final.data2, only.data)

# 158 - no censoring
only.data <- final.data1 %>% filter(Deer == "158")
only.data <- only.data %>% drop_na(Deer)
final.data2 <- rbind(final.data2, only.data)

# 160 - one relocation censored
only.data <- final.data1 %>% filter(Deer == "160")
only.data <- only.data %>% filter(UTMN < 4379500)
final.data2 <- rbind(final.data2, only.data)

# 161 - no censoring
only.data <- final.data1 %>% filter(Deer == "161")
final.data2 <- rbind(final.data2, only.data)

# 164 - no censoring
only.data <- final.data1 %>% filter(Deer == "164")
final.data2 <- rbind(final.data2, only.data)

# 166 - no censoring
only.data <- final.data1 %>% filter(Deer == "166")
final.data2 <- rbind(final.data2, only.data)

# 167 - no censoring
only.data <- final.data1 %>% filter(Deer == "167")
final.data2 <- rbind(final.data2, only.data)

# 168 - one relocation censored
only.data <- final.data1 %>% filter(Deer == "168")
only.data <- only.data %>% filter(UTME > 356200)
final.data2 <- rbind(final.data2, only.data)

# 169 - no censoring
only.data <- final.data1 %>% filter(Deer == "169")
final.data2 <- rbind(final.data2, only.data)

# 172 - no censoring
only.data <- final.data1 %>% filter(Deer == "172")
final.data2 <- rbind(final.data2, only.data)

# 173 - no censoring
only.data <- final.data1 %>% filter(Deer == "173")
final.data2 <- rbind(final.data2, only.data)

# 175 - no censoring
only.data <- final.data1 %>% filter(Deer == "175")
final.data2 <- rbind(final.data2, only.data)

# 176 - no censoring
only.data <- final.data1 %>% filter(Deer == "176")
final.data2 <- rbind(final.data2, only.data)

# 177 - no censoring
only.data <- final.data1 %>% filter(Deer == "177")
final.data2 <- rbind(final.data2, only.data)

# 178 - no censoring
only.data <- final.data1 %>% filter(Deer == "178")
final.data2 <- rbind(final.data2, only.data)

# 180 - no censoring
only.data <- final.data1 %>% filter(Deer == "180")
final.data2 <- rbind(final.data2, only.data)

# 181 - no censoring
only.data <- final.data1 %>% filter(Deer == "181")
final.data2 <- rbind(final.data2, only.data)

# 182 - no censoring
only.data <- final.data1 %>% filter(Deer == "182")
final.data2 <- rbind(final.data2, only.data)

# 192 - no censoring
only.data <- final.data1 %>% filter(Deer == "192")
final.data2 <- rbind(final.data2, only.data)

# 193 - several relocations censored
only.data <- final.data1 %>% filter(Deer == "193")
only.data <- only.data %>% filter(UTME < 358000)
final.data2 <- rbind(final.data2, only.data)

# 194 - no censoring
only.data <- final.data1 %>% filter(Deer == "194")
final.data2 <- rbind(final.data2, only.data)

# 195 - no censoring
only.data <- final.data1 %>% filter(Deer == "195")
final.data2 <- rbind(final.data2, only.data)

# 201 - no censoring
only.data <- final.data1 %>% filter(Deer == "201")
final.data2 <- rbind(final.data2, only.data)

# 202 - no censoring
only.data <- final.data1 %>% filter(Deer == "202")
final.data2 <- rbind(final.data2, only.data)

# 206 - no censoring
only.data <- final.data1 %>% filter(Deer == "206")
final.data2 <- rbind(final.data2, only.data)

# 207 - no censoring
only.data <- final.data1 %>% filter(Deer == "207")
final.data2 <- rbind(final.data2, only.data)

# 208 - no censoring
only.data <- final.data1 %>% filter(Deer == "208")
final.data2 <- rbind(final.data2, only.data)

# 209 - no censoring
only.data <- final.data1 %>% filter(Deer == "209")
final.data2 <- rbind(final.data2, only.data)

# 212 - no censoring
only.data <- final.data1 %>% filter(Deer == "212")
final.data2 <- rbind(final.data2, only.data)

# 213 - no censoring
only.data <- final.data1 %>% filter(Deer == "213")
final.data2 <- rbind(final.data2, only.data)

# 220 - some relocations censored
only.data <- final.data1 %>% filter(Deer == "220")
only.data <- only.data  %>% filter(UTME < 354000)
final.data2 <- rbind(final.data2, only.data)

# 221 - final relocation censored
only.data <- final.data1 %>% filter(Deer == "221")
only.data <- only.data[-5633,]
final.data2 <- rbind(final.data2, only.data)

# 222 - no censoring
only.data <- final.data1 %>% filter(Deer == "222")
final.data2 <- rbind(final.data2, only.data)

# 223 - no censoring
only.data <- final.data1 %>% filter(Deer == "223")
final.data2 <- rbind(final.data2, only.data)

# 226 - final relocation censored
only.data <- final.data1 %>% filter(Deer == "226")
only.data <- only.data[-6032,]
final.data2 <- rbind(final.data2, only.data)

# 227 - final relocations censored
only.data <- final.data1 %>% filter(Deer == "227")
only.data <- only.data[-c(6044:6047),]
final.data2 <- rbind(final.data2, only.data)

# 231 - no censoring
only.data <- final.data1 %>% filter(Deer == "231")
final.data2 <- rbind(final.data2, only.data)

# 233 - one relocation censored
only.data <- final.data1 %>% filter(Deer == "233")
only.data <- only.data  %>% filter(UTMN > 4378000)
final.data2 <- rbind(final.data2, only.data)

# 234 - several relocations censored
only.data <- final.data1 %>% filter(Deer == "234")
only.data <- only.data  %>% filter(UTME < 356000)
final.data2 <- rbind(final.data2, only.data)

# 236 - final relocation censored
only.data <- final.data1 %>% filter(Deer == "236")
only.data <- only.data[-5478,]
final.data2 <- rbind(final.data2, only.data)

# 237 - no censoring
only.data <- final.data1 %>% filter(Deer == "237")
final.data2 <- rbind(final.data2, only.data)

# 242 - final relocation censored
only.data <- final.data1 %>% filter(Deer == "242")
only.data <- only.data[-2466,]
final.data2 <- rbind(final.data2, only.data)

# 243 - no censoring
only.data <- final.data1 %>% filter(Deer == "243")
final.data2 <- rbind(final.data2, only.data)

# 244 - no censoring
only.data <- final.data1 %>% filter(Deer == "244")
final.data2 <- rbind(final.data2, only.data)

# 247 - no censoring
only.data <- final.data1 %>% filter(Deer == "247")
final.data2 <- rbind(final.data2, only.data)

# 250 - no censoring
only.data <- final.data1 %>% filter(Deer == "250")
final.data2 <- rbind(final.data2, only.data)

# 252 - no censoring
only.data <- final.data1 %>% filter(Deer == "252")
final.data2 <- rbind(final.data2, only.data)

# 255 - no censoring
only.data <- final.data1 %>% filter(Deer == "255")
final.data2 <- rbind(final.data2, only.data)

# 256 - no censoring
only.data <- final.data1 %>% filter(Deer == "256")
final.data2 <- rbind(final.data2, only.data)

# write to .csv
write.csv(final.data2, "deerdata_2.csv")

