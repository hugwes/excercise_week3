
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
ggplot(caro, aes(E, N, colour=static)) +
  geom_path(colour="black") +
  geom_point() +
  coord_fixed() +
  theme(legend.position="bottom") +
  theme_classic()
    
##########################################################
# Task 4: Segment-based Analysis


########################################################## 
# Task 5: Similarity Measures


########################################################## 
# Task 6: Calculate Similarity

  
  
  
  
  
