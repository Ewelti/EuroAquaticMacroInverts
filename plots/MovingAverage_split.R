##Set working directory
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/")

# attach data
MA <- read.csv("outputs/movingAve_split_YrEsts.csv")
MA$meanYr <- (MA$StartYear+ 4.5)
head(MA)

#plot for spp richness
SR <- subset(MA, Response == "SppRich")
SRs <- SR[ which(SR$meanYr > 1992), ]

SR_n <- subset(SRs, Lat == "north")
SR_s <- subset(SRs, Lat == "south")

tiff(filename = "MovingAverages_split.tiff", width = 11, height = 6, units = 'in', res = 600, compression = 'lzw')
par(mfrow=c(1,2),mar=c(4,4,0.4,0.4))

plot(SR_n$Estimate~SR_n$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-2.5,4))
title(ylab="Slope of species richness", line=2.4,cex.lab=1.3)
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70")
box(lwd=2)
points(SR_n$Estimate~c(SR_n$meanYr-0.25),pch="-", col="dodgerblue",cex=3)
points(SR_n$Estimate~c(SR_n$meanYr-0.25),type="l",lwd=2,col="dodgerblue")
arrows(c(SR_n$meanYr-0.25), SR_n$Q2.5, c(SR_n$meanYr-0.25), SR_n$Q97.5, length=0.05, angle=90, code=0,lwd=1,col="dodgerblue")
arrows(c(SR_n$meanYr-0.25), SR_n$Q5, c(SR_n$meanYr-0.25), SR_n$Q95, length=0.05, angle=90, code=0,lwd=3,col="dodgerblue")
arrows(c(SR_n$meanYr-0.25), SR_n$Q10, c(SR_n$meanYr-0.25), SR_n$Q90, length=0.05, angle=90, code=0,lwd=5,col="dodgerblue")

points(SR_s$Estimate~c(SR_s$meanYr+0.25),pch="-", col="brown",cex=3)
points(SR_s$Estimate~c(SR_s$meanYr+0.25),type="l",lwd=2,col="brown")
arrows(c(SR_s$meanYr+0.25), SR_s$Q2.5, c(SR_s$meanYr+0.25), SR_s$Q97.5, length=0.05, angle=90, code=0,lwd=1,col="brown")
arrows(c(SR_s$meanYr+0.25), SR_s$Q5, c(SR_s$meanYr+0.25), SR_s$Q95, length=0.05, angle=90, code=0,lwd=3,col="brown")
arrows(c(SR_s$meanYr+0.25), SR_s$Q10, c(SR_s$meanYr+0.25), SR_s$Q90, length=0.05, angle=90, code=0,lwd=5,col="brown")

#plot for abundance
AB <- subset(MA, Response == "abundance")
ABs <- AB[ which(AB$meanYr > 1992), ]

AB_n <- subset(ABs, Lat == "north")
AB_s <- subset(ABs, Lat == "south")

plot(AB_n$Estimate~AB_n$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-0.3,0.3))
title(ylab="Slope of abundance", line=2.4,cex.lab=1.3)
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70")
box(lwd=2)
points(AB_n$Estimate~c(AB_n$meanYr-0.25),pch="-", col="dodgerblue",cex=3)
points(AB_n$Estimate~c(AB_n$meanYr-0.25),type="l",lwd=2,col="dodgerblue")
arrows(c(AB_n$meanYr-0.25), AB_n$Q2.5, c(AB_n$meanYr-0.25), AB_n$Q97.5, length=0.05, angle=90, code=0,lwd=1,col="dodgerblue")
arrows(c(AB_n$meanYr-0.25), AB_n$Q5, c(AB_n$meanYr-0.25), AB_n$Q95, length=0.05, angle=90, code=0,lwd=3,col="dodgerblue")
arrows(c(AB_n$meanYr-0.25), AB_n$Q10, c(AB_n$meanYr-0.25), AB_n$Q90, length=0.05, angle=90, code=0,lwd=5,col="dodgerblue")

points(AB_s$Estimate~c(AB_s$meanYr+0.25),pch="-", col="brown",cex=3)
points(AB_s$Estimate~c(AB_s$meanYr+0.25),type="l",lwd=2,col="brown")
arrows(c(AB_s$meanYr+0.25), AB_s$Q2.5, c(AB_s$meanYr+0.25), AB_s$Q97.5, length=0.05, angle=90, code=0,lwd=1,col="brown")
arrows(c(AB_s$meanYr+0.25), AB_s$Q5, c(AB_s$meanYr+0.25), AB_s$Q95, length=0.05, angle=90, code=0,lwd=3,col="brown")
arrows(c(AB_s$meanYr+0.25), AB_s$Q10, c(AB_s$meanYr+0.25), AB_s$Q90, length=0.05, angle=90, code=0,lwd=5,col="brown")

legend("topright", legend=c("North","South"), pch=22, cex=2, col=c("dodgerblue","brown"),pt.bg=c("dodgerblue","brown"),bty="n")

dev.off()
