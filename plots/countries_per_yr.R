##Set working directory
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/")

# attach data
sites <- read.table("plots/sites_per_yr.csv", header=T, sep=",", row.names = 1, check.names = F)

#tiff(filename = "plots/Countries_yr.tiff", width = 8, height = 6.5, units = 'in', res = 600, compression = 'lzw')

par(mgp=c(3,0,-1),mar=c(3,5,1,1)+0.1)

colors <- c("darkblue", "red")
barplot(as.matrix(sites), las=2, col= colors)
legend(x=5,y=1300, rownames(sites), fill = colors, bty = "n")


dev.off()
##


