##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

# attach data
MA <- read.csv("outputs/SppLevelMovAve_YrEsts.csv")
MA$meanYr <- (MA$StartYear+ 4.5)
head(MA)
MA = subset(MA, select = -c(X) )
unique(MA$Response)
plot(MA$site_num~MA$StartYear)
abline(h=250)

#######################get data to get mean values of biodiversity metrics from
##attach data
sites <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T) # change file name according to the time series to be analyzed
#delete missing data rows
DATA2 <- sites[!is.na(sites$site_id_wMissing),]
#subset for sites with spp level data
sppD <- DATA2[which(DATA2$TaxonomicRes == "species"),]
sppD <- sppD[which(sppD$year >1993),]
length(unique(sppD$country))

tiff(filename = "plots/Fig3_movingWindow/SppLevel_MW.tiff", width = 8, height = 6.5, units = 'in', res = 600, compression = 'lzw')
par(mfrow=c(2,2),mar=c(2,4,0.2,0.2))

#plot for spp richness
SR <- subset(MA, Response == "spp_richness")
ave_SppRich <- mean(sppD$spp_richness)
st <-(SR[1:(nrow(SR)),3:12]/ave_SppRich)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1993 & sr$StartYear <2012), ] #full windows, min 200 sites

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
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
SRs <- sr[ which(sr$StartYear >1993 & sr$StartYear <2012), ] #full windows, min 200 sites

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
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

par(mar=c(4,4,0.2,0.2))
#plot for FRic
SR <- subset(MA, Response == "FRic")
st <-(SR[1:(nrow(SR)),3:12])*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1993 & sr$StartYear <2012), ] #full windows, min 200 sites

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
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
fr <-na.omit(sppD$FRed)
ave_FRed <- mean(fr)
st <-(SR[1:(nrow(SR)),3:12]/ave_FRed)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1993 & sr$StartYear <2012), ] #full windows, min 200 sites

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
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

