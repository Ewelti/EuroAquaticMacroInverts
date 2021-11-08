##Set working directory
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/")

# attach data
MA <- read.csv("outputs/movingAve_YrEsts.csv")
MA$meanYr <- (MA$StartYear+ 4.5)
head(MA)

tiff(filename = "MovingAverages.tiff", width = 11, height = 6, units = 'in', res = 600, compression = 'lzw')
par(mfrow=c(2,1),mar=c(4,4,0.4,0.4))

#plot for spp richness
SR <- subset(MA, Response == "spp_richness")
#select yrs with enough sites to be representative
#SRs <- SR[ which(SR$site_num >=50), ]
SRs <- SR[ which(SR$StartYear >=1989), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-0.5,1.3))#ylim=c(-0.8,2.55))
title(ylab="Slope of species richness", line=2.4,cex.lab=1.3)
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70")
box(lwd=2)
points(SRs$Estimate~SRs$meanYr,pch="-", col=1,cex=3)
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
arrows(SRs$meanYr, SRs$Q2.5, SRs$meanYr, SRs$Q97.5, length=0.05, angle=90, code=0,lwd=1)
arrows(SRs$meanYr, SRs$Q5, SRs$meanYr, SRs$Q95, length=0.05, angle=90, code=0,lwd=3)
arrows(SRs$meanYr, SRs$Q10, SRs$meanYr, SRs$Q90, length=0.05, angle=90, code=0,lwd=5)

#plot for abundance
SR <- subset(MA, Response == "abundance")
#select yrs with enough sites to be representative
#SRs <- SR[ which(SR$site_num >=50), ]
SRs <- SR[ which(SR$StartYear >=1989), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-0.06,0.06)) #ylim=c(-0.57,0.4))
title(ylab="Slope of abundance", line=2.7,cex.lab=1.3)
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70")
box(lwd=2)
points(SRs$Estimate~SRs$meanYr,pch="-", col=1,cex=3)
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
arrows(SRs$meanYr, SRs$Q2.5, SRs$meanYr, SRs$Q97.5, length=0.05, angle=90, code=0,lwd=1)
arrows(SRs$meanYr, SRs$Q5, SRs$meanYr, SRs$Q95, length=0.05, angle=90, code=0,lwd=3)
arrows(SRs$meanYr, SRs$Q10, SRs$meanYr, SRs$Q90, length=0.05, angle=90, code=0,lwd=5)

dev.off()


####
#####
par(mfrow=c(2,1),mar=c(2,4,0.4,0.4))

#plot for alien spp richness
SR <- subset(MA, Response == "alien_SppRich")
#select yrs with enough sites to be representative
SRs <- SR[ which(SR$site_num >=50), ]
#SRs <- SR[ which(SR$StartYear >=1989), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-0.04,0.1))
title(ylab="Slope of alien species richness", line=2.4,cex.lab=1.3)
title(xlab="", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70")
box(lwd=2)
points(SRs$Estimate~SRs$meanYr,pch="-", col=1,cex=3)
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
arrows(SRs$meanYr, SRs$Q2.5, SRs$meanYr, SRs$Q97.5, length=0.05, angle=90, code=0,lwd=1)
arrows(SRs$meanYr, SRs$Q5, SRs$meanYr, SRs$Q95, length=0.05, angle=90, code=0,lwd=3)
arrows(SRs$meanYr, SRs$Q10, SRs$meanYr, SRs$Q90, length=0.05, angle=90, code=0,lwd=5)

par(mar=c(4,4,0.4,0.4))
#plot for SppRich_nativeSpp
SR <- subset(MA, Response == "SppRich_nativeSpp")
#select yrs with enough sites to be representative
SRs <- SR[ which(SR$site_num >=50), ]
#SRs <- SR[ which(SR$StartYear >=1989), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-1.3,3.5)) #ylim=c(-0.57,0.4))
title(ylab="Slope of native spp richness", line=2.7,cex.lab=1.3)
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70")
box(lwd=2)
points(SRs$Estimate~SRs$meanYr,pch="-", col=1,cex=3)
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
arrows(SRs$meanYr, SRs$Q2.5, SRs$meanYr, SRs$Q97.5, length=0.05, angle=90, code=0,lwd=1)
arrows(SRs$meanYr, SRs$Q5, SRs$meanYr, SRs$Q95, length=0.05, angle=90, code=0,lwd=3)
arrows(SRs$meanYr, SRs$Q10, SRs$meanYr, SRs$Q90, length=0.05, angle=90, code=0,lwd=5)
