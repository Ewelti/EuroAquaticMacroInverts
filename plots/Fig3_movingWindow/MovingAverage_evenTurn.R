##Set working directory
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/")

# attach data
MA <- read.csv("outputs/movingAve_YrEsts.csv")
MA$meanYr <- (MA$StartYear+ 4.5)
head(MA)

tiff(filename = "plots/Fig3_movingWindow/MovingAverages_extra.tiff", width = 8, height = 6.5, units = 'in', res = 600, compression = 'lzw')
par(mfrow=c(2,2),mar=c(2,4,0.2,0.2))

#### Shannon's Evenness #####
SR <- subset(MA, Response == "E10")
st <-(SR[1:(nrow(SR)),3:12])*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
#SRs <- sr[ which(sr$site_num >=100), ]
SRs <- sr[ which(sr$StartYear >1989), ] # this is min 7 countries
#SRs <- sr[ which(sr$StartYear >=1980), ] #this is min 5 countries

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-2.5,2.5), xlim=c(1994.5,2018.5))#ylim=c(-0.8,2.55))
title(ylab=expression(paste("Shannon's evenness (% y"^"-1", ")")), line=2,cex.lab=1.4)
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

#plot for turnover
unique(MA$Response)
SR <- subset(MA, Response == "turnover")
ave_turnover <- 0.542933401
st <-(SR[1:(nrow(SR)),3:12]/ave_turnover)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1989), ] # this is min 7 countries

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-5.5,5.3), xlim=c(1994.5,2018.5)) #ylim=c(-0.57,0.4))
title(ylab=expression(paste("Turnover (% y"^"-1", ")")), line=2,cex.lab=1.4)
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

#### Func evenness #####
SR <- subset(MA, Response == "FEve")
ave_FEve <- 0.517910888
st <-(SR[1:(nrow(SR)),3:12]/ave_FEve)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >=1989), ] # this is min 7 countries

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-12,4), xlim=c(1994.5,2018.5)) #ylim=c(-0.57,0.4))
title(ylab=expression(paste("Functional evenness(% y"^"-1", ")")), line=2,cex.lab=1.4)
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

#### Func turnover #####
SR <- subset(MA, Response == "F_to")
st <-((exp(SR[1:(nrow(SR)),3:12])-1)*100)
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >=1989), ] # this is min 7 countries

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(-6,5), xlim=c(1994.5,2018.5)) #ylim=c(-0.57,0.4))
title(ylab=expression(paste("Functional turnover (% y"^"-1", ")")), line=2,cex.lab=1.4)
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
