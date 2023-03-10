##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

# attach data
response_stan <- readRDS("outputs/outputs_metaAnaylsis/stanTrends_site_level.rds")
head(response_stan)
unique(response_stan$Response)
length(unique(response_stan$site_id))

Ests <- read.csv("outputs/Yr_metaanaly_Ests.csv")

tiff(filename = "plots/Fig2_DensityPlots/SlopeDistributions_extra.tiff", width = 6, height = 6, units = 'in', res = 600, compression = 'lzw')

par(mar=c(4,0.4,0.4,0.4), mfrow=c(1,1))

#####################Extra metrics ##############################
#### Shannon's H #####
shannonsH <- subset(response_stan, Response == "shannonsH")
shannonsH <- shannonsH$estimate[!is.na(shannonsH$estimate)]
#average Shannon's H = 1.978232299
ave_shannonsH <- 1.978232299
percChange_perYr<-(shannonsH/ave_shannonsH)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
shannonsH_Est <- subset(Ests, Response == "shannonsH")
stand_shannonsH <- lapply(shannonsH_Est[,2:11],"*",100/ave_shannonsH)
yy <- (3.2/4*(b-a)+a)
points(x=stand_shannonsH$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_shannonsH$Q2.5, stand_shannonsH$Q2.5, stand_shannonsH$Q97.5, stand_shannonsH$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_shannonsH$Q5, stand_shannonsH$Q5, stand_shannonsH$Q95, stand_shannonsH$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_shannonsH$Q10, stand_shannonsH$Q10, stand_shannonsH$Q90, stand_shannonsH$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(4/4*(b-a)+a), legend=("a, Shannon's H"), bty="n", cex=1.3)
legend(x=7, y=(4/4*(b-a)+a), legend=expression(paste("+0.25% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### rarefied Spp Richness #####
SppRichRare <- subset(response_stan, Response == "spp_rich_rare")
SppRichRare <- SppRichRare$estimate[!is.na(SppRichRare$estimate)]
#average Rareified species richness= 19.2450259
ave_SppRichRare <- 19.2450259
percChange_perYr<-(SppRichRare/ave_SppRichRare)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
SppRichRare_Est <- subset(Ests, Response == "spp_richness_rarefied")
stand_SppRichRare <- lapply(SppRichRare_Est[,2:11],"*",100/ave_SppRichRare)
yy <- (2.2/4*(b-a)+a)
points(x=stand_SppRichRare$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_SppRichRare$Q2.5, stand_SppRichRare$Q2.5, stand_SppRichRare$Q97.5, stand_SppRichRare$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_SppRichRare$Q5, stand_SppRichRare$Q5, stand_SppRichRare$Q95, stand_SppRichRare$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_SppRichRare$Q10, stand_SppRichRare$Q10, stand_SppRichRare$Q90, stand_SppRichRare$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(3/4*(b-a)+a), legend=("b, Rarefied taxon richness"), bty="n", cex=1.3)
legend(x=7, y=(3/4*(b-a)+a), legend=expression(paste("+0.22% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Func diverg #####
FDiv <- subset(response_stan, Response == "FDiv")
FDiv <- FDiv$estimate[!is.na(FDiv$estimate)]
ave_FDiv <- 0.826257724
percChange_perYr<-(FDiv/(ave_FDiv^2))*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="Even", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
FDiv_Est <- subset(Ests, Response == "func_diverg")
stand_FDiv <- lapply(FDiv_Est[,2:11],"*",100/(ave_FDiv^2))
yy <- (1.2/4*(b-a)+a)
points(x=stand_FDiv$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FDiv$Q2.5, stand_FDiv$Q2.5, stand_FDiv$Q97.5, stand_FDiv$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FDiv$Q5, stand_FDiv$Q5, stand_FDiv$Q95, stand_FDiv$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FDiv$Q10, stand_FDiv$Q10, stand_FDiv$Q90, stand_FDiv$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(2/4*(b-a)+a), legend=("c, Func. divergence"), bty="n", cex=1.3)
legend(x=7, y=(2/4*(b-a)+a), legend=expression(paste("+0.09% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Rao's Q #####
RaoQ <- subset(response_stan, Response == "RaoQ")
RaoQ <- RaoQ$estimate[!is.na(RaoQ$estimate)]
#average Rao's Q = 37.7098828
ave_RaoQ <- 37.7098828
percChange_perYr<-(RaoQ/ave_RaoQ)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="Rao's Q", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
RaoQ_Est <- subset(Ests, Response == "RaoQ")
stand_RaoQ <- lapply(RaoQ_Est[,2:11],"*",100/ave_RaoQ)
yy <- (0.2/4*(b-a)+a)
points(x=stand_RaoQ$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_RaoQ$Q2.5, stand_RaoQ$Q2.5, stand_RaoQ$Q97.5, stand_RaoQ$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_RaoQ$Q5, stand_RaoQ$Q5, stand_RaoQ$Q95, stand_RaoQ$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_RaoQ$Q10, stand_RaoQ$Q10, stand_RaoQ$Q90, stand_RaoQ$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(1/4*(b-a)+a), legend=("d, Rao's Q"), bty="n", cex=1.3)
legend(x=7, y=(1/4*(b-a)+a), legend=expression(paste("+0.21% ",y^-1,, sep = "")), bty="n", cex=1.3)

##
box(lwd=2,col="white")
axis(1,lwd=2)
abline(v=0, lwd=1.5, lty=2)
##

dev.off()
########################################