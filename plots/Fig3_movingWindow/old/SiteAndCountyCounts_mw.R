##
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

## get mw data
ma <- read.csv("outputs/movingAve_YrEsts.csv")
ma$meanYr <- (ma$StartYear+ 4.5)
ma <- ma[which(ma$StartYear <2012), ]
head(ma)

## get htmw1 data
htmw1 <- read.csv("outputs/HighThresholdMovingAve1_YrEsts.csv")
htmw1$meanYr <- (htmw1$StartYear+ 4.5)
htmw1 <- htmw1[which(htmw1$StartYear <2010), ]
head(htmw1)

## get htmw2 data
htmw2 <- read.csv("outputs/HighThresholdMovingAve2_YrEsts.csv")
htmw2$meanYr <- (htmw2$StartYear+ 4.5)
htmw2 <- htmw2[which(htmw2$StartYear >1991 & htmw2$StartYear <2012), ]
head(htmw2)

## get spp level mw data
sppmw <- read.csv("outputs/SppLevelMovAve_YrEsts.csv")
sppmw$meanYr <- (sppmw$StartYear+ 4.5)
sppmw <- sppmw[which(sppmw$StartYear >1993 & sppmw$StartYear <2012), ]
head(sppmw)

#######################
tiff(filename = "plots/Fig3_movingWindow/MW_siteCount.tiff", width = 7, height = 6, units = 'in', res = 600, compression = 'lzw')

par(mfrow=c(2,2),mar=c(3,5,1,4))

ma$transco <- ma$country_num*70
max(ma$country_num)
plot(ma$site_num~ma$meanYr,type="n",ylab="", xlab="",las=1)
points(ma$site_num~ma$meanYr, type="l", col="black")
points(ma$site_num~ma$meanYr, pch=19)
points(ma$transco~ma$meanYr, type="l", col="cyan4")
points(ma$transco~ma$meanYr, pch=19, col="cyan4")
title(ylab="Number of sites", line=3, cex.lab=1.2)
title(xlab="Mean year of moving window", line=2, cex.lab=1.2)
legend("topleft", bty="n", legend=c("site count", "country count"), pch=19, col=c(1,"cyan4"))
legend("bottomright", bty="n", legend="a, Main analysis",cex=1.5)
marks <- c(8,10,12,14,16,18,20)
axis(4, at=c(560,700,840,980,1120,1260,1400),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of countries", line=2, las=3)
##

htmw1$transco <- htmw1$country_num*50
min(htmw1$country_num)
plot(htmw1$site_num~htmw1$meanYr,type="n",ylab="", xlab="",las=1)
points(htmw1$site_num~htmw1$meanYr, type="l", col="black")
points(htmw1$site_num~htmw1$meanYr, pch=19)
points(htmw1$transco~htmw1$meanYr, type="l", col="cyan4")
points(htmw1$transco~htmw1$meanYr, pch=19, col="cyan4")
title(ylab="Number of sites", line=3, cex.lab=1.2)
title(xlab="Mean year of moving window", line=2, cex.lab=1.2)
legend("bottomright", bty="n", legend="b, High threshold 1",cex=1.5)
marks <- c(11)
axis(4, at=c(550),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of countries", line=2, las=3)
##

htmw2$transco <- htmw2$country_num*35
max(htmw2$country_num)
plot(htmw2$site_num~htmw2$meanYr,type="n",ylab="", xlab="",las=1)
points(htmw2$site_num~htmw2$meanYr, type="l", col="black")
points(htmw2$site_num~htmw2$meanYr, pch=19)
points(htmw2$transco~htmw2$meanYr, type="l", col="cyan4")
points(htmw2$transco~htmw2$meanYr, pch=19, col="cyan4")
title(ylab="Number of sites", line=3, cex.lab=1.2)
title(xlab="Mean year of moving window", line=2, cex.lab=1.2)
legend("bottomright", bty="n", legend="c, High threshold 2",cex=1.5)
marks <- c(7,8)
axis(4, at=c(245,280),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of countries", line=2, las=3)
##

sppmw$transco <- sppmw$country_num*40
max(sppmw$country_num)
plot(sppmw$site_num~sppmw$meanYr,type="n",ylab="", xlab="",las=1)
points(sppmw$site_num~sppmw$meanYr, type="l", col="black")
points(sppmw$site_num~sppmw$meanYr, pch=19)
points(sppmw$transco~sppmw$meanYr, type="l", col="cyan4")
points(sppmw$transco~sppmw$meanYr, pch=19, col="cyan4")
title(ylab="Number of sites", line=3, cex.lab=1.2)
title(xlab="Mean year of moving window", line=2, cex.lab=1.2)
legend("bottomright", bty="n", legend="d, Spp level",cex=1.5)
marks <- c(6,8,10,12)
axis(4, at=c(240,320,400,480),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of countries", line=2, las=3)
##

dev.off()

##