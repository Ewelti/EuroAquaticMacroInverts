##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

# attach data
MA <- read.csv("outputs/outputs_movingWindow/movingAve_YrEsts.csv")
MA$meanYr <- (MA$StartYear+ 4.5)
head(MA)
MA = subset(MA, select = -c(X) )
unique(MA$Response)

metricdata <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv")
head(metricdata)

tiff(filename = "plots/Fig3_movingWindow/MovingAverages_extra.tif", width = 8, height = 6.5, units = 'in', res = 600, compression = 'lzw')
par(mfrow=c(2,2),mar=c(2,4,0.2,0.2))

#### Shannon's Evenness #####
SR <- subset(MA, Response == "E10")
st <-(10^(SR[1:(nrow(SR)),3:12])-1)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
#SRs <- sr[ which(sr$site_num >=100), ]
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ]
#SRs <- sr[ which(sr$StartYear >=1980), ] #this is min 5 countries

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
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
to<-as.numeric(replace(metricdata$turnover, metricdata$turnover=="<NA>", "NA"))
ave_turnover <- mean(to, na.rm=T)
st <-(SR[1:(nrow(SR)),3:12]/ave_turnover)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
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


par(mar=c(4,4,0.2,0.2))
#### Func evenness #####
SR <- subset(MA, Response == "FEve")
ave_FEve <- mean(metricdata$FEve, na.rm=T)
st <-(SR[1:(nrow(SR)),3:12]/ave_FEve)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
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
head(SR)

st <- SR[1:(nrow(SR)),3:12]

labs <- colnames(st)

backtrans <- function(x) {
perc <- ((plogis(interc + x * 5)/plogis(interc + x * -5))^(1/(10))-1)*100
print(perc)
}

for(i in 1:ncol(st)) {      
  st[,i] <- backtrans(st[,i])
}
colnames(st) <- labs

sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
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
