##Set working directory
setwd("C:/Users/Ellen/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

# attach data
sites <- read.table("plots/Fig1/sites_per_yr.csv", header=T, sep=",", row.names = 1, check.names = F)

tiff(filename = "plots/Fig1/Countries_yr.tiff", width =5, height = 8, units = 'in', res = 600, compression = 'lzw')

par(mgp=c(3,0,-1),mar=c(2,4,0.4,0.2)+0.1)

colors <- c("black","lightslateblue","lightskyblue", "darkblue","darkseagreen2", "forestgreen",
"brown1","darkorange","red3","peachpuff", "saddlebrown")

barplot(as.matrix(sites), las=2, col= colors, border=colors, xaxt='n')

place <- c(0.6,5.392,10.185,14.977,19.769,24.562,29.354,34.146,38.938,43.731,48.523,53.315,58.108,62.9)
marks <- c(1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)
axis(side = 1, at=place, labels=marks,las=2,tick=F,line=F)
legend(x=2,y=1400, rev(rownames(sites)), fill = rev(colors), border = rev(colors), bty = "n",cex=1.7)
title(ylab="Number of sites", line=2.2,cex.lab=1.5)

##
dev.off()
##


