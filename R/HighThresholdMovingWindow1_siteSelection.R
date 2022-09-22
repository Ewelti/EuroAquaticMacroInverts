setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")
#
rm(list=ls())

library(lubridate)

##attach data
sites <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T) # change file name according to the time series to be analyzed
#delete missing data rows
DATA2 <- sites[!is.na(sites$site_id_wMissing),]
sort(unique(DATA2$year))

#make turnover numeric
DATA2$turnover <- as.numeric(DATA2$turnover)

#count number of sites per year
sitecount <- aggregate(site_id ~ year, data = DATA2, FUN = length)
sitecount

#subset by year
sites_later1 <- subset(DATA2,year >= 2000)
sites_later <- subset(sites_later1,year <= 2018)
head(sites_later)

#count number of sampling years per site
yearcount <- aggregate(year ~ site_id, data = sites_later, FUN = length)
head(yearcount)

#subset for sites with more years
siteslater_long <- yearcount[which(yearcount$year >= 15),]
head(siteslater_long)
nrow(siteslater_long)

write.csv(siteslater_long, "sitesForHTMW1.csv")
