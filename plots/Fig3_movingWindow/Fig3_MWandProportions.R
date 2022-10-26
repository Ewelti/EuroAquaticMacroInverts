##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

#add pcks
library(tidyverse)
library(janitor)

# attach data
MA <- read.csv("outputs/movingAve_YrEsts.csv")
MA$meanYr <- (MA$StartYear+ 4.5)
head(MA)
MA = subset(MA, select = -c(X) )

#############################################

movav1 <- read.csv("outputs/MovingWindowProportions.csv")
##################################################

tiff(filename = "plots/Fig3_movingWindow/MovingAveragesProportions.tiff", width = 11, height = 6.5, units = 'in', res = 600, compression = 'lzw')

##layout
layout(mat = matrix(c(1:15), 
                        nrow = 3, 
                        ncol = 5),
       heights = c(2,2,0.4),
       widths = c(0.5,2,2,2,2))

par(mar=c(1.8,2,0.4,1))
#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Site % y"^"-1")), line=-2,cex.lab=1.6)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab="Proportion sites with", line=-0.3,cex.lab=1.4)
title(ylab="positive trends", line=-1.4,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

#plot for spp richness
SR <- subset(MA, Response == "spp_richness")
ave_SppRich <- 27.28712314
st <-(SR[1:(nrow(SR)),3:12]/ave_SppRich)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
#SRs <- sr[ which(sr$site_num >=100), ]
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ] # this is min 7 countries
#SRs <- sr[ which(sr$StartYear >=1980), ] #this is min 5 countries

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-1.4,3.7), xlim=c(1994.5,2015.5))
title(ylab=expression(paste("% y"^"-1")), line=2.4,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)#pch="-",
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("topright", bty="n", legend="a, Taxon richness",cex=1.5)

count_sign_mw <- movav1[which(movav1$response=='spp_richness'),]
plot(0,0, ylim=c(0.4,0.7),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
legend("topright", bty="n", legend="e, Taxon richness",cex=1.5)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-2,cex.lab=1.3)

#plot for abundance
SR <- subset(MA, Response == "abundance")
st <-(SR[1:(nrow(SR)),3:12]*100)
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-5.3,3.6), xlim=c(1994.5,2015.5))
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)#pch="-",
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("topright", bty="n", legend="b, Abundance",cex=1.5)

count_sign_mw <- movav1[which(movav1$response=='abundance'),]
plot(0,0, ylim=c(0.4,0.7),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
title(xlab="Mean year of moving window", line=2.5,cex.lab=1.1)
legend("topright", bty="n", legend="f, Abundance",cex=1.5)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-2,cex.lab=1.3)

#plot for FRic
SR <- subset(MA, Response == "FRic")
st <-(SR[1:(nrow(SR)),3:12])*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-2,4), xlim=c(1994.5,2015.5))
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("topright", bty="n", legend="c, Functional richness",cex=1.5)

count_sign_mw <- movav1[which(movav1$response=='FRic'),]
plot(0,0, ylim=c(0.4,0.7),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
title(xlab="Mean year of moving window", line=2.5,cex.lab=1.1)
legend("topright", bty="n", legend="g, Functional richness",cex=1.5)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-2,cex.lab=1.3)

#plot for FRed
SR <- subset(MA, Response == "FRed")
ave_FRed <- 0.291145846
st <-(SR[1:(nrow(SR)),3:12]/ave_FRed)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-1,2), xlim=c(1994.5,2015.5))
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("topright", bty="n", legend="d, Functional redundancy",cex=1.5)

count_sign_mw <- movav1[which(movav1$response=='FRed'),]
plot(0,0, ylim=c(0.4,0.7),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
title(xlab="Mean year of moving window", line=2.5,cex.lab=1.1)
legend("topright", bty="n", legend="h, Functional redundancy",cex=1.5)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-2,cex.lab=1.3)

##

dev.off()
###############################################################################


