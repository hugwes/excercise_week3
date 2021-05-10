
### EXCERCISE 3 ###

##########################################################
# Preparations
# https://github.com/hugwes/excercise_week3 (github-link)

##########################################################
# Load in Libraries
library(readr)        # to import tabular data (e.g. csv)
library(dplyr)        # to manipulate (tabular) data
library(ggplot2)      # to visualize data
library(sf)           # to handle spatial vector data
library(terra)        # to handle raster data
library(lubridate)    # to handle dates and times
library(zoo)          # moving window functions

##########################################################
# Import Data
caro <- read_delim("caro60.csv",",")

##########################################################
# Task 1: Segmentation

# Specify a Temporal Window for in which to measure Euclidean Distances
# The sampling interval is 1 minute! 
# We take a Temporal Window of 6 minutes!
# That would mean including 6 fixes! (-1,-2,-3,+1,+2,+3) 

# Measure the distance from every point to every other point within this Temporal Window
caro <- caro %>%
  mutate(
    n_minus1 = sqrt((lag(E,1)-E)^2+(lag(N,1)-N)^2),     # distance to pos -1 minutes
    n_minus2 = sqrt((lag(E,2)-E)^2+(lag(N,2)-N)^2),     # distance to pos -2 minutes
    n_minus3 = sqrt((lag(E,3)-E)^2+(lag(N,3)-N)^2),     # distance to pos -3 minutes
    n_plus1 = sqrt((E-lead(E,1))^2+(N-lead(N,1))^2),    # distance to pos +1 mintues
    n_plus2 = sqrt((E-lead(E,2))^2+(N-lead(N,2))^2),    # distance to pos +2 minutes
    n_plus3 = sqrt((E-lead(E,3))^2+(N-lead(N,3))^2))    # distance to pos +3 minutes

caro <- caro %>%
  rowwise() %>%
  mutate(stepMean = mean(c(n_minus1, n_minus2, n_minus3, n_plus1, n_plus2, n_plus3))) %>%
  ungroup()

##########################################################
# Task 2: Specify and Apply Threshold d
# Exploring the values
summary(caro$stepMean)
hist(caro$stepMean)
boxplot(caro$stepMean)

# Useing the mean of all stepMean values as threshold
# Generating a new column with the name static (stops = TRUE, moves = FALSE)
caro <- caro %>% 
  ungroup() %>%
  mutate(static = stepMean < mean(stepMean, na.rm = TRUE))

##########################################################
# Task 3: Visualize Segmented Trajectories
ggplot(caro, aes(E,N,colour=static)) +
  geom_path(colour="black") +
  geom_point() +
  coord_fixed() +
  theme(legend.position="bottom") +
  theme_classic()
    
##########################################################
# Task 4: Segment-based Analysis
# Create a function to assign unique IDs based on the column static
rle_id <- function(vec){
  x <- rle(vec)$lengths
  as.factor(rep(seq_along(x), times=x))}

# Add a new coloum with the unique ID's 
caro <- caro %>%
  mutate(segment_id = rle_id(static))

# We are just interested in the moves and not in the stops
caro_filter_1 <- caro %>%
  filter(!static)

# Plot 1: Moving segments coloured by the segment_id (all segments)
ggplot(caro_filter_1, aes(E,N,colour=segment_id)) +
  geom_path() +
  geom_point() +
  coord_fixed() +
  labs(title="all segments") +
  theme_classic()

# Plot 2: Moving segments coloured by the segment_id (long segments)
caro_filter_2 <- caro %>%
  filter(!static,
         segment_id != "8",
         segment_id != "10",
         segment_id != "12")       

ggplot(caro_filter_2, aes(E,N,colour=segment_id)) +
  geom_path() +
  geom_point() +
  coord_fixed() +
  labs(title="long segments") +
  theme_classic()

########################################################## 
# Task 5: Similarity Measures
# Import Data
pedestrian <- read_delim("pedestrian.csv",",")

ggplot(pedestrian, aes(E,N,colour=TrajID)) +
  geom_path() +
  geom_point() +
  coord_fixed() +
  labs(title="Visual comparison of the 6 trajectories") +
  facet_wrap(~TrajID) + 
  theme_classic()

########################################################## 
# Task 6: Calculate Similarity
# Get to know the library SimilarityMeasures
library(SimilarityMeasures)
help(package = "SimilarityMeasures")

# Filter by Trajectory & Convert them into Matrix
ped1 <- pedestrian %>% 
  filter(TrajID == 1)
ped1 <- matrix(data=c(as.numeric(ped1$E),as.numeric(ped1$N)), ncol = 2)

ped2 <- pedestrian %>%
  filter(TrajID == 2)
ped2 <- matrix(data=c(as.numeric(ped2$E),as.numeric(ped2$N)), ncol = 2)

ped3 <- pedestrian %>%
  filter(TrajID == 3) 
ped3 <- matrix(data=c(as.numeric(ped3$E), as.numeric(ped3$N)), ncol=2)

ped4 <- pedestrian %>%
  filter(TrajID == 4) 
ped4 <- matrix(data=c(as.numeric(ped4$E), as.numeric(ped4$N)), ncol=2)

ped5 <- pedestrian %>%
  filter(TrajID == 5) 
ped5 <- matrix(data=c(as.numeric(ped5$E), as.numeric(ped5$N)), ncol=2)

ped6 <- pedestrian %>%
  filter(TrajID == 6) 
ped6 <- matrix(data=c(as.numeric(ped6$E), as.numeric(ped6$N)), ncol=2)

# Dynamic Time Warping (DTW)
DTW1_2 <- DTW(ped1, ped2)
DTW1_3 <- DTW(ped1, ped3)
DTW1_4 <- DTW(ped1, ped4)
DTW1_5 <- DTW(ped1, ped5)
DTW1_6 <- DTW(ped1, ped6)

# Edit Distance (EditDist)
ED1_2 <- EditDist(ped1, ped2)
ED1_3 <- EditDist(ped1, ped3)
ED1_4 <- EditDist(ped1, ped4)
ED1_5 <- EditDist(ped1, ped5)
ED1_6 <- EditDist(ped1, ped6)

# Frechet Calculation (Frechnet)
Frechet1_2 <- Frechet(ped1, ped2)
Frechet1_3 <- Frechet(ped1, ped3)
Frechet1_4 <- Frechet(ped1, ped4)
Frechet1_5 <- Frechet(ped1, ped5)
Frechet1_6 <- Frechet(ped1, ped6)

# Longest Common Subsequence (LCSS) --> takes ages!!!
# Frechet1_2 <- LCSS(ped1, ped2)
# Frechet1_3 <- LCSS(ped1, ped3)
# Frechet1_4 <- LCSS(ped1, ped4)
# Frechet1_5 <- LCSS(ped1, ped5)
# Frechet1_6 <- LCSS(ped1, ped6)

# Add those Values back to the Dataframe (pedestrian)
pedestrian <- pedestrian %>%
  mutate(DTW = case_when(
    TrajID == "1" ~ 0,
    TrajID == "2" ~ DTW1_2,
    TrajID == "3" ~ DTW1_3,
    TrajID == "4" ~ DTW1_4,
    TrajID == "5" ~ DTW1_5,
    TrajID == "6" ~ DTW1_6)) %>%
  mutate(EditDist = case_when(
    TrajID == "1" ~ 0,
    TrajID == "2" ~ ED1_2,
    TrajID == "3" ~ ED1_3,
    TrajID == "4" ~ ED1_4,
    TrajID == "5" ~ ED1_5,
    TrajID == "6" ~ ED1_6)) %>%
  mutate(Frechet = case_when(
    TrajID == "1" ~ 0,
    TrajID == "2" ~ Frechet1_2,
    TrajID == "3" ~ Frechet1_3,
    TrajID == "4" ~ Frechet1_4,
    TrajID == "5" ~ Frechet1_5,
    TrajID == "6" ~ Frechet1_6))

# Plot 1: Dynamic Time Warping (DTW)
ggplot(pedestrian, aes(TrajID, DTW, fill=TrajID)) +
  geom_col() + 
  labs(title = "DTW")
  
# Plot 2: Edit Distance
ggplot(pedestrian, aes(TrajID, EditDist, fill=TrajID)) +
  geom_col() + 
  labs(title = "EditDist")

# Plot 3: Frechet
ggplot(pedestrian, aes(TrajID, Frechet, fill=TrajID)) +
  geom_col() + 
  labs(title = "Frechet")



