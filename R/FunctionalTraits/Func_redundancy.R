#set wd
setwd("C:/Users/ewelti/Desktop/aquatic_data/TraitCalculation")

#libraries
#install.packages("FD")
library(ade4)
library(adiv)
library(vegan)
library(FD)

## load data
# Comm data
head(comm <- read.csv("mform_abundances.csv", h = T, sep = ",", row.names = 1, stringsAsFactors = FALSE)) # site by species matrix
# Trait data
head(traits <- read.csv("Traits_fuzzy_BiolOnly_85percCoverage.csv", h = T, sep = ",", row.names = 1, stringsAsFactors = FALSE)) # species by trait matrix

# Check if row names and column names match
identical(row.names(traits), colnames(comm))

# Create distance matrix based on traits
(dis <- dist(traits, method = "euclidean")) # Create distance matrix
(dis <- as.dist(decostand(dis, "range"))) # Values in dis must be bounded between 0 and 1.
(gowdis <- gowdis(traits)) # test gowdis because I think it is used for the other FD metrics (FRic, FEve, FDis, etc.)

# Calculate BAT package FD metrics at the community level
(func_div <- uniqueness(comm, dis, tol = 1e-08, abundance = TRUE)) # creates 3 dataframs: kbar, V, and Red
str(func_div)
# Isolate redundancy metrics
redund <- func_div$red # these are calculated at the community level
head(redund)
# Save redundancy dataframe
write.csv(redund, "Comm_redund.csv")
