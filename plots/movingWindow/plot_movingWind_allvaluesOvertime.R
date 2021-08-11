##Set working directory
setwd("C:/Users/ewelti/Desktop/aquatic_data/movingWindow")
#setwd(dir = "C:/Users/Ellen/Desktop/aquatic_data/movingWindow")


# attach data
go <- read.csv("AllMovingWindow.csv", header=T) # change file name according to the time series to be analyzed
attach(go)
head(go)

par(mfrow=c(2,3),mar=c(2,4,0.4,0.4))
goo <- go[!is.na(sppRich),]
plot(goo$sppRich~goo$meanYear,ylab="Slope of species richness +/- SE", cex=1.5, xlab="", type="n", las=1, ylim=c(-0.2,0.7))
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70") #abline(h=0,lty=2)
points(goo$sppRich~goo$meanYear,pch=19, col=1,cex=1.5)
points(goo$sppRich~goo$meanYear,type="l",lwd=2)
arrows(meanYear, sppRich-sppRich_se, meanYear, sppRich+sppRich_se, length=0.05, angle=90, code=3,lwd=2)
box()

goo <- go[!is.na(abun),]
plot(goo$abun~goo$meanYear,ylab="Slope of abundance +/- SE", cex=1.5, xlab="", type="n", las=1, ylim=c(-0.03,0.07))
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70") #abline(h=0,lty=2)
points(goo$abun~goo$meanYear,pch=19, col=1,cex=1.5)
points(goo$abun~goo$meanYear,type="l",lwd=2)
arrows(meanYear, abun-abun_se, meanYear, abun+abun_se, length=0.05, angle=90, code=3,lwd=2)
box()

goo <- go[!is.na(ShanH),]
plot(goo$ShanH~goo$meanYear,ylab="Slope of Shannon's diversity +/- SE", cex=1.5, xlab="", type="n", las=1, ylim=c(-0.01,0.04))
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70") #abline(h=0,lty=2)
points(goo$ShanH~goo$meanYear,pch=19, col=1,cex=1.5)
points(goo$ShanH~goo$meanYear,type="l",lwd=2)
arrows(meanYear, ShanH-ShanH_se, meanYear, ShanH+ShanH_se, length=0.05, angle=90, code=3,lwd=2)
box()

goo <- go[!is.na(sppRichRare),]
plot(goo$sppRichRare~goo$meanYear,ylab="Slope of rarefied species richnes +/- SE", cex=1.5, xlab="", type="n", las=1, ylim=c(-0.2,0.4))
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70") #abline(h=0,lty=2)
points(goo$sppRichRare~goo$meanYear,pch=19, col=1,cex=1.5)
points(goo$sppRichRare~goo$meanYear,type="l",lwd=2)
arrows(meanYear, sppRichRare-sppRichRare_se, meanYear, sppRichRare+sppRichRare_se, length=0.05, angle=90, code=3,lwd=2)
box()

goo <- go[!is.na(S_PIE),]
plot(goo$S_PIE~goo$meanYear,ylab="Slope of S_PIE +/- SE", cex=1.5, xlab="", type="n", las=1, ylim=c(-0.01,0.02))
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70") #abline(h=0,lty=2)
points(goo$S_PIE~goo$meanYear,pch=19, col=1,cex=1.5)
points(goo$S_PIE~goo$meanYear,type="l",lwd=2)
arrows(meanYear, S_PIE-S_PIE_se, meanYear, S_PIE+S_PIE_se, length=0.05, angle=90, code=3,lwd=2)
box()

goo <- go[!is.na(turnover),]
plot(goo$turnover~goo$meanYear,ylab="Slope of turnover +/- SE", cex=1.5, xlab="", type="n", las=1, ylim=c(-0.025,0.015))
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey70",border = "grey70") #abline(h=0,lty=2)
points(goo$turnover~goo$meanYear,pch=19, col=1,cex=1.5)
points(goo$turnover~goo$meanYear,type="l",lwd=2)
arrows(meanYear, turnover-turnover_se, meanYear, turnover+turnover_se, length=0.05, angle=90, code=3,lwd=2)
box()