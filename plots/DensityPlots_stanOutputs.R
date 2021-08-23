##Set working directory
#setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/")

# attach data
response_stan <- readRDS("outputs/stanTrends_site_level.rds")
head(response_stan)
unique(response_stan$Response)

Ests <- read.csv("outputs/Yr_metaanaly_Ests.csv")

#attach data from second run
response_stan_logged <- readRDS("outputs/stanTrends_site_level_logged.rds")
Ests_logged <- read.csv("outputs/Yr_metaanaly_Ests_logged.csv")

par(mar=c(4,0.4,0.4,0.4), mfrow=c(1,3))

#####################Taxonomic metrics ##############################
#### Spp Richness #####
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
#axis(2, at=0, labels="SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
SppRich_Est <- subset(Ests, Response == "spp_richness")
stand_SppRich <- lapply(SppRich_Est[,2:11],"*",100/ave_SppRich)
yy <- (5/6*(b-a)+a)
points(x=stand_SppRich$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_SppRich$Q2.5, stand_SppRich$Q2.5, stand_SppRich$Q97.5, stand_SppRich$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_SppRich$Q5, stand_SppRich$Q5, stand_SppRich$Q95, stand_SppRich$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(6/6*(b-a)+a), legend=("Species richness"), bty="n", cex=1.3)

#### rarefied Spp Richness #####
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
#axis(2, at=0, labels="Sn", las=1,cex.axis=1.3)
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
polygon(x=c(stand_SppRichRare$Q2.5, stand_SppRichRare$Q2.5, stand_SppRichRare$Q97.5, stand_SppRichRare$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_SppRichRare$Q5, stand_SppRichRare$Q5, stand_SppRichRare$Q95, stand_SppRichRare$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(5/6*(b-a)+a), legend=("Rarefied spp richness"), bty="n", cex=1.3)

#### Abundance #####
abund <- subset(response_stan, Response == "abundance")
abund <- abund$estimate[!is.na(abund$estimate)]
percChange_perYr <- abund*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="Abun", las=1,cex.axis=1.3)
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
polygon(x=c(stand_abund$Q2.5, stand_abund$Q2.5, stand_abund$Q97.5, stand_abund$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_abund$Q5, stand_abund$Q5, stand_abund$Q95, stand_abund$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(4/6*(b-a)+a), legend=("Abundance"), bty="n", cex=1.3)

#### Shannon's H #####
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
#axis(2, at=0, labels="ShanH", las=1,cex.axis=1.3)
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
polygon(x=c(stand_shannonsH$Q2.5, stand_shannonsH$Q2.5, stand_shannonsH$Q97.5, stand_shannonsH$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_shannonsH$Q5, stand_shannonsH$Q5, stand_shannonsH$Q95, stand_shannonsH$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(3/6*(b-a)+a), legend=("Shannon's H"), bty="n", cex=1.3)

#### Shannon's Evenness #####
E10 <- subset(response_stan_logged, Response == "E10")
E10 <- E10$estimate[!is.na(E10$estimate)]
#average Shannon's Evenness = 0.085019802
percChange_perYr<-(E10*100)
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*5
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="Even", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
E10_Est <- subset(Ests_logged, Response == "E10")
stand_E10 <- lapply(E10_Est[,2:11],"*",100)
yy <- (1/6*(b-a)+a)
points(x=stand_E10$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_E10$Q2.5, stand_E10$Q2.5, stand_E10$Q97.5, stand_E10$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_E10$Q5, stand_E10$Q5, stand_E10$Q95, stand_E10$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(2/6*(b-a)+a), legend=("Shannon's evenness"), bty="n", cex=1.3)

#### Turnover #####
turnover <- subset(response_stan_logged, Response == "turnover")
turnover <- turnover$estimate[!is.na(turnover$estimate)]
head(turnover)
#average turnover = 0.542933401
ave_turnover <- 0.542933401
percChange_perYr<-(turnover/ave_turnover)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*6
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="TurnO", las=1,cex.axis=1.3)
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
polygon(x=c(stand_turnover$Q2.5, stand_turnover$Q2.5, stand_turnover$Q97.5, stand_turnover$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_turnover$Q5, stand_turnover$Q5, stand_turnover$Q95, stand_turnover$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(1/6*(b-a)+a), legend=("Turnover"), bty="n", cex=1.3)

##
box(lwd=2)
abline(v=0, lwd=1, lty=2)
##

#####################Functional metrics ##############################
#### Functional Richness #####
FRic <- subset(response_stan_logged, Response == "spp_richness")
FRic <- FRic$estimate[!is.na(FRic$estimate)]
#average functional richness= 1372.61497
percChange_perYr<-(FRic*100)
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-4
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n')
title(xlab="Annual % change", line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="FRic", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
FRic_Est <- subset(Ests_logged, Response == "func_rich")
stand_FRic <- lapply(FRic_Est[,2:11],"*", 100)
yy <- (4/5*(b-a)+a)
points(x=stand_FRic$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_FRic$Q2.5, stand_FRic$Q2.5, stand_FRic$Q97.5, stand_FRic$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_FRic$Q5, stand_FRic$Q5, stand_FRic$Q95, stand_FRic$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(5/5*(b-a)+a), legend=("Func. richness"), bty="n", cex=1.3)

#### Rao's Q #####
RaoQ <- subset(response_stan, Response == "RaoQ")
RaoQ <- RaoQ$estimate[!is.na(RaoQ$estimate)]
#average Rao's Q = 37.7098828
ave_RaoQ <- 37.7098828
percChange_perYr<-(RaoQ/ave_RaoQ)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="Rao's Q", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
RaoQ_Est <- subset(Ests, Response == "RaoQ")
stand_RaoQ <- lapply(RaoQ_Est[,2:11],"*",100/ave_RaoQ)
yy <- (3/5*(b-a)+a)
points(x=stand_RaoQ$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_RaoQ$Q2.5, stand_RaoQ$Q2.5, stand_RaoQ$Q97.5, stand_RaoQ$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_RaoQ$Q5, stand_RaoQ$Q5, stand_RaoQ$Q95, stand_RaoQ$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(4/5*(b-a)+a), legend=("Rao's Q"), bty="n", cex=1.3)

#### Func diverg #####
FDiv <- subset(response_stan_logged, Response == "FDiv")
FDiv <- FDiv$estimate[!is.na(FDiv$estimate)]
#average FDiv = 0.826257724
ave_FDiv <- 0.826257724
percChange_perYr<-(FDiv/(ave_FDiv^2))*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="FDiv", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
FDiv_Est <- subset(Ests_logged, Response == "func_diverg")
stand_FDiv <- lapply(FDiv_Est[,2:11],"*",100/(ave_FDiv^2))
yy <- (2/5*(b-a)+a)
points(x=stand_FDiv$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_FDiv$Q2.5, stand_FDiv$Q2.5, stand_FDiv$Q97.5, stand_FDiv$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_FDiv$Q5, stand_FDiv$Q5, stand_FDiv$Q95, stand_FDiv$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(3/5*(b-a)+a), legend=("Func. divergence"), bty="n", cex=1.3)

#### Func evenness #####
FEve <- subset(response_stan, Response == "FEve")
FEve <- FEve$estimate[!is.na(FEve$estimate)]
#average Func even = 0.517910888
ave_FEve <- 0.517910888
percChange_perYr<-(FEve/ave_FEve)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-1
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="FEve", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
FEve_Est <- subset(Ests, Response == "func_even")
stand_FEve <- lapply(FEve_Est[,2:11],"*",100/ave_FEve)
yy <- (1/5*(b-a)+a)
points(x=stand_FEve$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_FEve$Q2.5, stand_FEve$Q2.5, stand_FEve$Q97.5, stand_FEve$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_FEve$Q5, stand_FEve$Q5, stand_FEve$Q95, stand_FEve$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(2/5*(b-a)+a), legend=("Func. evenness"), bty="n", cex=1.3)

#### Func turnover #####
F_to <- subset(response_stan_logged, Response == "F_to")
F_to <- F_to$estimate[!is.na(F_to$estimate)]
#average Func turnover = 0.086496081
#ave_F_to <- 0.086496081
percChange_perYr<-(F_to*100)
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*5
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="FTurno", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
F_to_Est <- subset(Ests_logged, Response == "func_turnover")
stand_F_to <- lapply(F_to_Est[,2:11],"*",100)
yy <- (0/5*(b-a)+a)
points(x=stand_F_to$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_F_to$Q2.5, stand_F_to$Q2.5, stand_F_to$Q97.5, stand_F_to$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_F_to$Q5, stand_F_to$Q5, stand_F_to$Q95, stand_F_to$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(1/5*(b-a)+a), legend=("Func. turnover"), bty="n", cex=1.3)

##
box(lwd=2)
abline(v=0, lwd=1, lty=2)
##

#####################Subsets ##############################
#### Alien Spp Richness #####
alien_SppRich <- subset(response_stan, Response == "alien_SppRich")
alien_SppRich <- alien_SppRich$estimate[!is.na(alien_SppRich$estimate)]
#average alien richness= 1.420037807
ave_alien_SppRich <- 1.420037807
percChange_perYr<-(alien_SppRich/ave_alien_SppRich)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-7
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n')
title(xlab="Annual % change", line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="Alien SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
alien_SppRich_Est <- subset(Ests, Response == "alien_SppRich")
stand_alien_SppRich <- lapply(alien_SppRich_Est[,2:11],"*", 100/ave_alien_SppRich)
yy <- (7/8*(b-a)+a)
points(x=stand_alien_SppRich$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_alien_SppRich$Q2.5, stand_alien_SppRich$Q2.5, stand_alien_SppRich$Q97.5, stand_alien_SppRich$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_alien_SppRich$Q5, stand_alien_SppRich$Q5, stand_alien_SppRich$Q95, stand_alien_SppRich$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(8/8*(b-a)+a), legend=("Alien spp richness"), bty="n", cex=1.3)

#### alien abund #####
alien_Abund <- subset(response_stan, Response == "alien_Abund")
alien_Abund <- alien_Abund$estimate[!is.na(alien_Abund$estimate)]
percChange_perYr<- alien_Abund*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-6
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="Alien ab", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
alien_Abund_Est <- subset(Ests, Response == "alien_abund")
stand_alien_Abund <- lapply(alien_Abund_Est[,2:11],"*",100)
yy <- (6/8*(b-a)+a)
points(x=stand_alien_Abund$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_alien_Abund$Q2.5, stand_alien_Abund$Q2.5, stand_alien_Abund$Q97.5, stand_alien_Abund$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_alien_Abund$Q5, stand_alien_Abund$Q5, stand_alien_Abund$Q95, stand_alien_Abund$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(7/8*(b-a)+a), legend=("Alien abundance"), bty="n", cex=1.3)

#### native spp rich #####
SppRich_nativeSpp <- subset(response_stan, Response == "SppRich_nativeSpp")
SppRich_nativeSpp <- SppRich_nativeSpp$estimate[!is.na(SppRich_nativeSpp$estimate)]
#average native spp rich = 27.52824362
ave_SppRich_nativeSpp <- 27.52824362
percChange_perYr<-(SppRich_nativeSpp/ave_SppRich_nativeSpp)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-5
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="Native SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
SppRich_nativeSpp_Est <- subset(Ests, Response == "native_SppRich")
stand_SppRich_nativeSpp <- lapply(SppRich_nativeSpp_Est[,2:11],"*",100/ave_SppRich_nativeSpp)
yy <- (5/8*(b-a)+a)
points(x=stand_SppRich_nativeSpp$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_SppRich_nativeSpp$Q2.5, stand_SppRich_nativeSpp$Q2.5, stand_SppRich_nativeSpp$Q97.5, stand_SppRich_nativeSpp$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_SppRich_nativeSpp$Q5, stand_SppRich_nativeSpp$Q5, stand_SppRich_nativeSpp$Q95, stand_SppRich_nativeSpp$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(6/8*(b-a)+a), legend=("Native spp richness"), bty="n", cex=1.3)

#### native abund #####
abund_nativeSpp <- subset(response_stan, Response == "abund_nativeSpp")
abund_nativeSpp <- abund_nativeSpp$estimate[!is.na(abund_nativeSpp$estimate)]
percChange_perYr<- abund_nativeSpp*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-4
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="Native ab", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
abund_nativeSpp_Est <- subset(Ests, Response == "native_abund")
stand_abund_nativeSpp <- lapply(abund_nativeSpp_Est[,2:11],"*",100)
yy <- (4/8*(b-a)+a)
points(x=stand_abund_nativeSpp$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_abund_nativeSpp$Q2.5, stand_abund_nativeSpp$Q2.5, stand_abund_nativeSpp$Q97.5, stand_abund_nativeSpp$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_abund_nativeSpp$Q5, stand_abund_nativeSpp$Q5, stand_abund_nativeSpp$Q95, stand_abund_nativeSpp$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(5/8*(b-a)+a), legend=("Native abundance"), bty="n", cex=1.3)

#### EPT spp rich #####
EPT_SppRich <- subset(response_stan, Response == "EPT_SppRich")
EPT_SppRich <- EPT_SppRich$estimate[!is.na(EPT_SppRich$estimate)]
#average EPT spp rich = 9.279548523
ave_EPT_SppRich <- 9.279548523
percChange_perYr<-(EPT_SppRich/ave_EPT_SppRich)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*5
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="EPT SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
EPT_SppRich_Est <- subset(Ests, Response == "EPT_SppRich")
stand_EPT_SppRich <- lapply(EPT_SppRich_Est[,2:11],"*",100/ave_EPT_SppRich)
yy <- (3/8*(b-a)+a)
points(x=stand_EPT_SppRich$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_EPT_SppRich$Q2.5, stand_EPT_SppRich$Q2.5, stand_EPT_SppRich$Q97.5, stand_EPT_SppRich$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_EPT_SppRich$Q5, stand_EPT_SppRich$Q5, stand_EPT_SppRich$Q95, stand_EPT_SppRich$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(4/8*(b-a)+a), legend=("EPT spp richness"), bty="n", cex=1.3)

#### EPT abund #####
EPT_Abund <- subset(response_stan, Response == "EPT_Abund")
EPT_Abund <- EPT_Abund$estimate[!is.na(EPT_Abund$estimate)]
percChange_perYr<- EPT_Abund*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*6
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="EPT ab", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
EPT_Abund_Est <- subset(Ests, Response == "EPT_abund")
stand_EPT_Abund <- lapply(EPT_Abund_Est[,2:11],"*",100)
yy <- (2/8*(b-a)+a)
points(x=stand_EPT_Abund$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_EPT_Abund$Q2.5, stand_EPT_Abund$Q2.5, stand_EPT_Abund$Q97.5, stand_EPT_Abund$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_EPT_Abund$Q5, stand_EPT_Abund$Q5, stand_EPT_Abund$Q95, stand_EPT_Abund$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(3/8*(b-a)+a), legend=("EPT abundance"), bty="n", cex=1.3)

#### insect spp rich #####
insect_SppRich <- subset(response_stan, Response == "insect_SppRich")
insect_SppRich <- insect_SppRich$estimate[!is.na(insect_SppRich$estimate)]
#average insect spp rich = 19.19656517
ave_insect_SppRich <- 19.19656517
percChange_perYr<-(insect_SppRich/ave_insect_SppRich)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-1
b <- (max(d$y)+(max(d$y)/10))*7
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="Insect SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
insect_SppRich_Est <- subset(Ests, Response == "insect_SppRich")
stand_insect_SppRich <- lapply(insect_SppRich_Est[,2:11],"*",100/ave_insect_SppRich)
yy <- (1/8*(b-a)+a)
points(x=stand_insect_SppRich$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_insect_SppRich$Q2.5, stand_insect_SppRich$Q2.5, stand_insect_SppRich$Q97.5, stand_insect_SppRich$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_insect_SppRich$Q5, stand_insect_SppRich$Q5, stand_insect_SppRich$Q95, stand_insect_SppRich$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(2/8*(b-a)+a), legend=("Insect spp richness"), bty="n", cex=1.3)

#### insect abund #####
insect_Abund <- subset(response_stan, Response == "insect_Abund")
insect_Abund <- insect_Abund$estimate[!is.na(insect_Abund$estimate)]
percChange_perYr<- insect_Abund*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*8
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="Insect ab", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
insect_Abund_Est <- subset(Ests, Response == "insect_abund")
stand_insect_Abund <- lapply(insect_Abund_Est[,2:11],"*",100)
yy <- (0/8*(b-a)+a)
points(x=stand_insect_Abund$Estimate, y=yy, lwd=2,pch="l",cex=1.5)
polygon(x=c(stand_insect_Abund$Q2.5, stand_insect_Abund$Q2.5, stand_insect_Abund$Q97.5, stand_insect_Abund$Q97.5),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 2)
polygon(x=c(stand_insect_Abund$Q5, stand_insect_Abund$Q5, stand_insect_Abund$Q95, stand_insect_Abund$Q95),
        y=c(yy,yy,yy,yy),
        col = 1,border = 1,lwd = 3)
legend(x=-12, y=(1/8*(b-a)+a), legend=("Insect abundance"), bty="n", cex=1.3)
box(lwd=2)
##
abline(v=0, lwd=1, lty=2)
##
########################################