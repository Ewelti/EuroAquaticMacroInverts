

# attach data
response_stan <- readRDS("outputs/stanTrends_site_level.rds")
head(response_stan)
unique(response_stan$Response)


Ests <- read.csv("outputs/Yr_metaanaly_Ests.csv")

par(mar=c(4,5,0.4,0.4))
##
SppRich <- subset(response_stan, Response == "spp_richness")
SppRich <- SppRich$estimate[!is.na(SppRich$estimate)]
#average species richness= 27.28712314 species
ave_SppRich <- 27.28712314
percChange_perYr<-(SppRich/ave_SppRich)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-5
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n')
title(xlab="Annual % change", line=2.4,cex.lab=1.3)
axis(2, at=0, labels="Spp R", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
SppRich_Est <- subset(Ests, Response == "spp_richness")
stand_SppRich <- lapply(SppRich_Est[,2:11],"*",100/ave_sppRich)
yy <- (5/6*(b-a)+a)
points(x=stand_SppRich$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_SppRich$Q1, stand_SppRich$Q1, stand_SppRich$Q99, stand_SppRich$Q99),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_SppRich$Q5, stand_SppRich$Q5, stand_SppRich$Q95, stand_SppRich$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 4)

##
SppRichRare <- subset(response_stan, Response == "spp_rich_rare")
SppRichRare <- SppRichRare$estimate[!is.na(SppRichRare$estimate)]
#average Rareified species richness= 19.2450259
ave_SppRichRare <- 19.2450259
percChange_perYr<-(SppRichRare/ave_SppRichRare)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-4
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
axis(2, at=0, labels="Sn", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
SppRichRare_Est <- subset(Ests, Response == "spp_richness_rarefied")
stand_SppRichRare <- lapply(SppRichRare_Est[,2:11],"*",100/ave_SppRichRare)
yy <- (4/6*(b-a)+a)
points(x=stand_SppRichRare$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_SppRichRare$Q1, stand_SppRichRare$Q1, stand_SppRichRare$Q99, stand_SppRichRare$Q99),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_SppRichRare$Q5, stand_SppRichRare$Q5, stand_SppRichRare$Q95, stand_SppRichRare$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 4)

##
abund <- subset(response_stan, Response == "abundance")
abund <- abund$estimate[!is.na(abund$estimate)]
ab<-Abun_Est[!is.na(Abun_Est)]
percChange_perYr <- abund*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
axis(2, at=0, labels="Abun", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
abund_Est <- subset(Ests, Response == "abundance")
stand_abund <- lapply(abund_Est[,2:11],"*",100)
yy <- (3/6*(b-a)+a)
points(x=stand_abund$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_abund$Q1, stand_abund$Q1, stand_abund$Q99, stand_abund$Q99),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_abund$Q5, stand_abund$Q5, stand_abund$Q95, stand_abund$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 4)

##
shannonsH <- subset(response_stan, Response == "shannonsH")
shannonsH <- shannonsH$estimate[!is.na(shannonsH$estimate)]
#average Shannon's H = 1.978232299
ave_shannonsH <- 1.978232299
percChange_perYr<-(shannonsH/ave_shannonsH)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
axis(2, at=0, labels="ShanH", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
shannonsH_Est <- subset(Ests, Response == "shannonsH")
stand_shannonsH <- lapply(shannonsH_Est[,2:11],"*",100/ave_shannonsH)
yy <- (2/6*(b-a)+a)
points(x=stand_shannonsH$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_shannonsH$Q1, stand_shannonsH$Q1, stand_shannonsH$Q99, stand_shannonsH$Q99),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_shannonsH$Q5, stand_shannonsH$Q5, stand_shannonsH$Q95, stand_shannonsH$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 4)

#### Shannon's Evenness #####
E10 <- subset(response_stan, Response == "E10")
E10 <- E10$estimate[!is.na(E10$estimate)]
#average Shannon's Evenness = 0.085019802
ave_E10 <- 0.085019802
percChange_perYr<-(E10/ave_E10)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-1
b <- (max(d$y)+(max(d$y)/10))*5
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
axis(2, at=0, labels="Even", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
E10_Est <- subset(Ests, Response == "E10")
stand_E10 <- lapply(E10_Est[,2:11],"*",100/ave_E10)
yy <- (1/6*(b-a)+a)
points(x=stand_E10$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_E10$Q1, stand_E10$Q1, stand_E10$Q99, stand_E10$Q99),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_E10$Q5, stand_E10$Q5, stand_E10$Q95, stand_E10$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 4)

#### Turnover #####
turnover <- subset(response_stan, Response == "turnover")
turnover <- turnover$estimate[!is.na(turnover$estimate)]
#average turnover = 0.542933401
ave_turnover <- 0.542933401
percChange_perYr<-(turnover/ave_turnover)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*6
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
axis(2, at=0, labels="TurnO", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
turnover_Est <- subset(Ests, Response == "turnover")
stand_turnover <- lapply(turnover_Est[,2:11],"*",100/ave_turnover)
yy <- (0/6*(b-a)+a)
points(x=stand_turnover$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_turnover$Q1, stand_turnover$Q1, stand_turnover$Q99, stand_turnover$Q99),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_turnover$Q5, stand_turnover$Q5, stand_turnover$Q95, stand_turnover$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 4)

##
abline(v=0, lwd=1, lty=2)
##