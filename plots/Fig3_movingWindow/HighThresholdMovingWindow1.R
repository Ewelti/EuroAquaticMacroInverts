##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

# attach data
MA <- read.csv("outputs/HighThresholdMovingAve1_YrEsts.csv")
MA$meanYr <- (MA$StartYear+ 4.5)
head(MA)
MA = subset(MA, select = -c(X) )
unique(MA$Response)
plot(MA$site_num~MA$StartYear)
abline(h=500)

tiff(filename = "plots/Fig3_movingWindow/HTMW1.tiff", width = 8, height = 6.5, units = 'in', res = 600, compression = 'lzw')
par(mfrow=c(2,2),mar=c(2,4,0.2,0.2))

#plot for spp richness
SR <- subset(MA, Response == "spp_richness")
ave_SppRich <- 27.28712314
st <-(SR[1:(nrow(SR)),3:12]/ave_SppRich)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1999 & sr$StartYear <2010), ] # this is min 500 sites and within selected years of high number of sites

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-1.5,1.5), xlim=c(2004.5,2013.5))#ylim=c(-0.8,2.55))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
#title(xlab="Mean year of moving window", line=2,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
abline(h=0,lty=2)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
#box(lwd=2)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)#pch="-",
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("topright", bty="n", legend="a",cex=1.5)

#plot for abundance
SR <- subset(MA, Response == "abundance")
st <-(SR[1:(nrow(SR)),3:12]*100)
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1999 & sr$StartYear <2010), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-2.7,2.5), xlim=c(2004.5,2013.5)) #ylim=c(-0.57,0.4))
title(ylab=expression(paste("Abundance (% y"^"-1", ")")), line=2,cex.lab=1.4)
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
abline(h=0,lty=2)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)#pch="-",
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("topright", bty="n", legend="b",cex=1.5)

#plot for FRic
SR <- subset(MA, Response == "FRic")
st <-(SR[1:(nrow(SR)),3:12])*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1999 & sr$StartYear <2010), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-3.6,2.5), xlim=c(2004.5,2013.5)) #ylim=c(-0.57,0.4))
title(ylab=expression(paste("Functional richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
abline(h=0,lty=2)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)#pch="-",
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("topright", bty="n", legend="c",cex=1.5)

#plot for FRed
SR <- subset(MA, Response == "FRed")
ave_FRed <- 0.291145846
st <-(SR[1:(nrow(SR)),3:12]/ave_FRed)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1999 & sr$StartYear <2010), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-1.3,1), xlim=c(2004.5,2013.5)) #ylim=c(-0.57,0.4))
title(ylab=expression(paste("Functional redundancy (% y"^"-1", ")")), line=2,cex.lab=1.4)
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
abline(h=0,lty=2)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)#pch="-",
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("topright", bty="n", legend="d",cex=1.5)

dev.off()


####
#####################
###############################################
################################################

