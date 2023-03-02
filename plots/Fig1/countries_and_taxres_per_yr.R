##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

# attach data
sites <- read.table("plots/Fig1/sites_per_yr.csv", header=T, sep=",", row.names = 1, check.names = F)

tiff(filename = "plots/Fig1/Countries_yr.tiff", width =5, height = 8, units = 'in', res = 600, compression = 'lzw')

par(mgp=c(3,0,-1),mar=c(2,4,0.4,0.2)+0.1)

colors <- c("black","lightslateblue","lightskyblue", "darkblue","darkseagreen2", "forestgreen",
"brown1","darkorange","red3","peachpuff", "saddlebrown")

barplot(as.matrix(sites), las=2, col= colors, border=colors, xaxt='n',space=0)

place <- c(0.6,5.392,10.185,14.977,19.769,24.562,29.354,34.146,38.938,43.731,48.523,53.315,58.108,62.9)
marks <- c(1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)
axis(side = 1, at=place, labels=marks,las=2,tick=F,line=F)
legend(x=2,y=1400, rev(rownames(sites)), fill = rev(colors), border = rev(colors), bty = "n",cex=1.7)
title(ylab="Number of sites", line=2.2,cex.lab=1.5)

##
dev.off()
##
d1 <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T) # change file name according to the time series to be analyzed
allYrs <- d1[!is.na(d1$site_id_wMissing),]
head(allYrs)

d1_or <- allYrs[order(allYrs$year),]
DT <- data.table(d1_or)

# number of sites per taxonomic resolution per year
TaxResCount <- DT[, .(site_num = length(unique(site_id))), by = list(year,TaxonomicRes)]
head(TaxResCount)

TaxResCount <- TaxResCount[order(TaxResCount$year),]


###########sites of a given taxonomic Res count/year
library(tidyr)
trw <- spread(TaxResCount, year, site_num)
trw <- as.data.frame(trw)
trw [is.na(trw )] <- 0

# # bargraph # #
tiff(filename = "plots/TaxaResPerYr.tiff", width =8, height = 8, units = 'in', res = 600, compression = 'lzw')

par(mgp=c(3,0,-1),mar=c(4,4,0.4,0.2)+0.1)

colors <- c("black","lightskyblue","lightslateblue")

barplot(as.matrix(trw), las=2, col= colors, border=colors, xaxt='n')

place <- c(0.6,5.392,10.185,14.977,19.769,24.562,29.354,34.146,38.938,43.731,48.523,53.315,58.108,62.9)
marks <- c(1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)
axis(side = 1, at=place, labels=marks,las=2,tick=F,line=F)
legend(x=2,y=1200, rev(trw[,1]), fill = rev(colors), border = rev(colors), bty = "n",cex=1.7)
title(ylab="Number of sites", line=2.2,cex.lab=1.5)

##
dev.off()
##

