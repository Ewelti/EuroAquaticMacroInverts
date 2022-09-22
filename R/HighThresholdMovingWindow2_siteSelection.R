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
sites_later <- subset(DATA2,year > 1989)
head(sites_later)

#count number of sampling years per site
yearcount <- aggregate(year ~ site_id, data = sites_later, FUN = length)
head(yearcount)

#subset for sites with more years
siteslater_long <- yearcount[which(yearcount$year >= 20),]
head(siteslater_long)
nrow(siteslater_long)

write.csv(siteslater_long, "yearsSampled_since1990.csv")

#subset full dataset for sites with at least yrs of data in 1990 or later
allYrs <- subset(sites_later, site_id %in% siteslater_long$site_id)
head(allYrs)

#see how many sites have each number of years of data 8-30 years from 1990 onward
HTMWc <- NULL
for(i in 8:30){
	siteslater_long <- yearcount[which(yearcount$year >= i),]
	numsites <- nrow(siteslater_long)
	d.i <- data.frame(years = i, site_count = numsites)
	HTMWc <- rbind(HTMWc, d.i)
}
HTMWc