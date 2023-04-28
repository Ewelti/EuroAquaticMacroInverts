##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

# attach data
response_stan <- readRDS("outputs/outputs_metaAnalysis/stanTrends_site_level.rds")
head(response_stan)
unique(response_stan$Response)
length(unique(response_stan$site_id))

metricdata <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv")
head(metricdata)

Ests <- read.csv("outputs/outputs_metaAnalysis/Yr_metaanaly_Ests.csv")

tiff(filename = "plots/Fig2_DensityPlots/SlopeDistributions_Subsets.tif", width = 12, height = 6, units = 'in', res = 600, compression = 'lzw')

par(mar=c(4,0.4,0.4,0.4), mfrow=c(1,2))

#####################Taxonomic metrics ##############################
#####################Subsets ##############################
#### Alien Spp Richness #####
alien_SppRich <- subset(response_stan, Response == "alien_SppRich")
alien_SppRich <- alien_SppRich$estimate[!is.na(alien_SppRich$estimate)]
md <- metricdata[ which(metricdata$alien_detected=='yes'), ]
ave_alien_SppRich <- mean(md$alien_SppRich, na.rm=T)
percChange_perYr<-(alien_SppRich/ave_alien_SppRich)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n', yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="Alien SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
alien_SppRich_Est <- subset(Ests, Response == "alien_SppRich_unweighted")
stand_alien_SppRich <- lapply(alien_SppRich_Est[,2:11],"*", (100/ave_alien_SppRich))
yy <- (3.2/4*(b-a)+a)
points(x=stand_alien_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_alien_SppRich$Q2.5, stand_alien_SppRich$Q2.5, stand_alien_SppRich$Q97.5, stand_alien_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_alien_SppRich$Q5, stand_alien_SppRich$Q5, stand_alien_SppRich$Q95, stand_alien_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_alien_SppRich$Q10, stand_alien_SppRich$Q10, stand_alien_SppRich$Q90, stand_alien_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(4/4*(b-a)+a), legend=("a, Non-native richness"), bty="n", cex=1.3)
legend(x=7, y=(4/4*(b-a)+a), legend=expression(paste("+3.97% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### native spp rich #####
SppRich_nativeSpp <- subset(response_stan, Response == "SppRich_nativeSpp")
SppRich_nativeSpp <- SppRich_nativeSpp$estimate[!is.na(SppRich_nativeSpp$estimate)]
ave_SppRich_nativeSpp <- mean(metricdata$SppRich_nativeSpp, na.rm=T)
percChange_perYr<-(SppRich_nativeSpp/ave_SppRich_nativeSpp)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n', yaxs="i",xaxt='n')
#axis(2, at=0, labels="Native SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
SppRich_nativeSpp_Est <- subset(Ests, Response == "native_SppRich")
stand_SppRich_nativeSpp <- lapply(SppRich_nativeSpp_Est[,2:11],"*",(100/ave_SppRich_nativeSpp))
yy <- (2.2/4*(b-a)+a)
points(x=stand_SppRich_nativeSpp$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_SppRich_nativeSpp$Q2.5, stand_SppRich_nativeSpp$Q2.5, stand_SppRich_nativeSpp$Q97.5, stand_SppRich_nativeSpp$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_SppRich_nativeSpp$Q5, stand_SppRich_nativeSpp$Q5, stand_SppRich_nativeSpp$Q95, stand_SppRich_nativeSpp$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_SppRich_nativeSpp$Q10, stand_SppRich_nativeSpp$Q10, stand_SppRich_nativeSpp$Q90, stand_SppRich_nativeSpp$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(3/4*(b-a)+a), legend=("c, Native richness"), bty="n", cex=1.3)
legend(x=7, y=(3/4*(b-a)+a), legend=expression(paste("+0.64% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### EPT spp rich #####
EPT_SppRich <- subset(response_stan, Response == "EPT_SppRich")
EPT_SppRich <- EPT_SppRich$estimate[!is.na(EPT_SppRich$estimate)]
ave_EPT_SppRich <- mean(metricdata$EPT_SppRich, na.rm=T)
percChange_perYr<-(EPT_SppRich/ave_EPT_SppRich)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-1
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n', yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="EPT SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
EPT_SppRich_Est <- subset(Ests, Response == "EPT_SppRich_unweighted")
stand_EPT_SppRich <- lapply(EPT_SppRich_Est[,2:11],"*",100/ave_EPT_SppRich)
yy <- (1.2/4*(b-a)+a)
points(x=stand_EPT_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_EPT_SppRich$Q2.5, stand_EPT_SppRich$Q2.5, stand_EPT_SppRich$Q97.5, stand_EPT_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_EPT_SppRich$Q5, stand_EPT_SppRich$Q5, stand_EPT_SppRich$Q95, stand_EPT_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_EPT_SppRich$Q10, stand_EPT_SppRich$Q10, stand_EPT_SppRich$Q90, stand_EPT_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(2/4*(b-a)+a), legend=("e, EPT richness"), bty="n", cex=1.3)
legend(x=7, y=(2/4*(b-a)+a), legend=expression(paste("+0.45% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### insect spp rich #####
insect_SppRich <- subset(response_stan, Response == "insect_SppRich")
insect_SppRich <- insect_SppRich$estimate[!is.na(insect_SppRich$estimate)]
ave_insect_SppRich <- mean(metricdata$insect_SppRich, na.rm=T)
percChange_perYr<-(insect_SppRich/ave_insect_SppRich)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n', yaxs="i",xaxt='n')
#axis(2, at=0, labels="Insect SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
insect_SppRich_Est <- subset(Ests, Response == "insect_SppRich")
stand_insect_SppRich <- lapply(insect_SppRich_Est[,2:11],"*",100/ave_insect_SppRich)
yy <- (0.2/4*(b-a)+a)
points(x=stand_insect_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_insect_SppRich$Q2.5, stand_insect_SppRich$Q2.5, stand_insect_SppRich$Q97.5, stand_insect_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_insect_SppRich$Q5, stand_insect_SppRich$Q5, stand_insect_SppRich$Q95, stand_insect_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_insect_SppRich$Q10, stand_insect_SppRich$Q10, stand_insect_SppRich$Q90, stand_insect_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(1/4*(b-a)+a), legend=("g, Insect richness"), bty="n", cex=1.3)
legend(x=7, y=(1/4*(b-a)+a), legend=expression(paste("+0.71% ",y^-1,, sep = "")), bty="n", cex=1.3)
##
box(lwd=2,col="white")
axis(1,lwd=2)
abline(v=0, lwd=1.5, lty=2)

##
#### alien abund #####
alien_Abund <- subset(response_stan, Response == "alien_Abund")
alien_Abund <- alien_Abund$estimate[!is.na(alien_Abund$estimate)]
percChange_perYr<- (10^alien_Abund-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-3
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n', yaxs="i",xaxt='n')
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="Alien ab", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
alien_Abund_Est <- subset(Ests, Response == "alien_abund")
stand_alien_Abund <- lapply((10^alien_Abund_Est[,2:11]-1),"*",100)
yy <- (3.2/4*(b-a)+a)
points(x=stand_alien_Abund$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_alien_Abund$Q2.5, stand_alien_Abund$Q2.5, stand_alien_Abund$Q97.5, stand_alien_Abund$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_alien_Abund$Q5, stand_alien_Abund$Q5, stand_alien_Abund$Q95, stand_alien_Abund$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_alien_Abund$Q10, stand_alien_Abund$Q10, stand_alien_Abund$Q90, stand_alien_Abund$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(4/4*(b-a)+a), legend=("b, Non-native abundance"), bty="n", cex=1.3)
legend(x=7, y=(4/4*(b-a)+a), legend=expression(paste("+3.9% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### native abund #####
abund_nativeSpp <- subset(response_stan, Response == "abund_nativeSpp")
abund_nativeSpp <- abund_nativeSpp$estimate[!is.na(abund_nativeSpp$estimate)]
percChange_perYr<- (10^abund_nativeSpp-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n', yaxs="i",xaxt='n')
#axis(2, at=0, labels="Native ab", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
abund_nativeSpp_Est <- subset(Ests, Response == "native_abund")
stand_abund_nativeSpp <- lapply((10^abund_nativeSpp_Est[,2:11]-1),"*",100)
yy <- (2.2/4*(b-a)+a)
points(x=stand_abund_nativeSpp$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_abund_nativeSpp$Q2.5, stand_abund_nativeSpp$Q2.5, stand_abund_nativeSpp$Q97.5, stand_abund_nativeSpp$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_abund_nativeSpp$Q5, stand_abund_nativeSpp$Q5, stand_abund_nativeSpp$Q95, stand_abund_nativeSpp$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_abund_nativeSpp$Q10, stand_abund_nativeSpp$Q10, stand_abund_nativeSpp$Q90, stand_abund_nativeSpp$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(3/4*(b-a)+a), legend=("d, Native abundance"), bty="n", cex=1.3)
legend(x=7, y=(3/4*(b-a)+a), legend=expression(paste("+0.26% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### EPT abund #####
EPT_Abund <- subset(response_stan, Response == "EPT_Abund")
EPT_Abund <- EPT_Abund$estimate[!is.na(EPT_Abund$estimate)]
percChange_perYr<- (10^EPT_Abund-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-1
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n', yaxs="i",xaxt='n')
#axis(2, at=0, labels="EPT ab", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
EPT_Abund_Est <- subset(Ests, Response == "EPT_abund_unweighted")
stand_EPT_Abund <- lapply((10^EPT_Abund_Est[,2:11]-1),"*",100)
yy <- (1.2/4*(b-a)+a)
points(x=stand_EPT_Abund$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_EPT_Abund$Q2.5, stand_EPT_Abund$Q2.5, stand_EPT_Abund$Q97.5, stand_EPT_Abund$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_EPT_Abund$Q5, stand_EPT_Abund$Q5, stand_EPT_Abund$Q95, stand_EPT_Abund$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_EPT_Abund$Q10, stand_EPT_Abund$Q10, stand_EPT_Abund$Q90, stand_EPT_Abund$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(2/4*(b-a)+a), legend=("f, EPT abundance"), bty="n", cex=1.3)
legend(x=7, y=(2/4*(b-a)+a), legend=expression(paste("+2.38% ",y^-1,, sep = "")), bty="n", cex=1.3)


#### insect abund #####
insect_Abund <- subset(response_stan, Response == "insect_Abund")
insect_Abund <- insect_Abund$estimate[!is.na(insect_Abund$estimate)]
percChange_perYr<- (10^insect_Abund-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n', yaxs="i",xaxt='n')
#axis(2, at=0, labels="Insect ab", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
insect_Abund_Est <- subset(Ests, Response == "insect_abund")
stand_insect_Abund <- lapply((10^insect_Abund_Est[,2:11]-1),"*",100)
yy <- (0.2/4*(b-a)+a)
points(x=stand_insect_Abund$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_insect_Abund$Q2.5, stand_insect_Abund$Q2.5, stand_insect_Abund$Q97.5, stand_insect_Abund$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_insect_Abund$Q5, stand_insect_Abund$Q5, stand_insect_Abund$Q95, stand_insect_Abund$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_insect_Abund$Q10, stand_insect_Abund$Q10, stand_insect_Abund$Q90, stand_insect_Abund$Q90),
        y=c((yy-yy/4.5),(yy+yy/4.5),(yy+yy/4.5),(yy-yy/4.6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(1/4*(b-a)+a), legend=("h, Insect abundance"), bty="n", cex=1.3)
legend(x=7, y=(1/4*(b-a)+a), legend=expression(paste("+1.53% ",y^-1,, sep = "")), bty="n", cex=1.3)
##
box(lwd=2,col="white")
axis(1,lwd=2)
abline(v=0, lwd=1.5, lty=2)
##

dev.off()
########################################