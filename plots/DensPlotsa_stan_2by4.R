##Set working directory
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/")

# attach data
response_stan <- readRDS("outputs/stanTrends_site_level.rds")
head(response_stan)
unique(response_stan$Response)
length(unique(response_stan$site_id))

Ests <- read.csv("outputs/Yr_metaanaly_Ests.csv")

par(mar=c(4,0.4,0.4,0.4), mfrow=c(1,2))

#####################Taxonomic metrics ##############################
#### Abundance #####
abund <- subset(response_stan, Response == "abundance")
abund <- abund$estimate[!is.na(abund$estimate)]
percChange_perYr <- abund*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
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
yy <- (3.2/4*(b-a)+a)
points(x=stand_abund$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_abund$Q2.5, stand_abund$Q2.5, stand_abund$Q97.5, stand_abund$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_abund$Q5, stand_abund$Q5, stand_abund$Q95, stand_abund$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_abund$Q10, stand_abund$Q10, stand_abund$Q90, stand_abund$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(4/4*(b-a)+a), legend=("Abundance"), bty="n", cex=1.3)
legend(x=5, y=(4/4*(b-a)+a), legend=expression(paste("+1.13% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Spp Richness #####
SppRich <- subset(response_stan, Response == "spp_richness")
SppRich <- SppRich$estimate[!is.na(SppRich$estimate)]
#average species richness= 27.28712314 species
ave_SppRich <- 27.28712314
percChange_perYr<-(SppRich/ave_SppRich)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n')
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
yy <- (2.2/4*(b-a)+a)
points(x=stand_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_SppRich$Q2.5, stand_SppRich$Q2.5, stand_SppRich$Q97.5, stand_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_SppRich$Q5, stand_SppRich$Q5, stand_SppRich$Q95, stand_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_SppRich$Q10, stand_SppRich$Q10, stand_SppRich$Q90, stand_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd = 1)
legend(x=-18, y=(3/4*(b-a)+a), legend=("Species richness"), bty="n", cex=1.3)
legend(x=5, y=(3/4*(b-a)+a), legend=expression(paste("+0.94% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Shannon's Evenness #####
E10 <- subset(response_stan, Response == "E10")
E10 <- E10$estimate[!is.na(E10$estimate)]
#average Shannon's Evenness = 0.085019802
percChange_perYr<-(E10*100)
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="Even", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
E10_Est <- subset(Ests, Response == "E10")
stand_E10 <- lapply(E10_Est[,2:11],"*",100)
yy <- (1.2/4*(b-a)+a)
points(x=stand_E10$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_E10$Q2.5, stand_E10$Q2.5, stand_E10$Q97.5, stand_E10$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_E10$Q5, stand_E10$Q5, stand_E10$Q95, stand_E10$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_E10$Q10, stand_E10$Q10, stand_E10$Q90, stand_E10$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(2/4*(b-a)+a), legend=("Shannon's evenness"), bty="n", cex=1.3)
legend(x=5, y=(2/4*(b-a)+a), legend=expression(paste("-0.5% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Turnover #####
turnover <- subset(response_stan, Response == "turnover")
turnover <- turnover$estimate[!is.na(turnover$estimate)]
head(turnover)
#average turnover = 0.542933401
ave_turnover <- 0.542933401
percChange_perYr<-(turnover/ave_turnover)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
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
yy <- (0.2/4*(b-a)+a)
points(x=stand_turnover$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_turnover$Q2.5, stand_turnover$Q2.5, stand_turnover$Q97.5, stand_turnover$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_turnover$Q5, stand_turnover$Q5, stand_turnover$Q95, stand_turnover$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_turnover$Q10, stand_turnover$Q10, stand_turnover$Q90, stand_turnover$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(1/4*(b-a)+a), legend=("Turnover"), bty="n", cex=1.3)
legend(x=5, y=(1/4*(b-a)+a), legend=expression(paste("-0.3% ",y^-1,, sep = "")), bty="n", cex=1.3)

##
box(lwd=2)
abline(v=0, lwd=1, lty=2)
##

#####################Functional metrics ##############################
#### Functional Redundancy #####
FRed <- subset(response_stan, Response == "FRed")
FRed <- FRed$estimate[!is.na(FRed$estimate)]
#average Rao's Q = 0.291145846
ave_FRed <- 0.291145846
percChange_perYr<-(FRed/ave_FRed)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="Rao's Q", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
FRed_Est <- subset(Ests, Response == "FRed")
stand_FRed <- lapply(FRed_Est[,2:11],"*",100/ave_FRed)
yy <- (3.2/4*(b-a)+a)
points(x=stand_FRed$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FRed$Q2.5, stand_FRed$Q2.5, stand_FRed$Q97.5, stand_FRed$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FRed$Q5, stand_FRed$Q5, stand_FRed$Q95, stand_FRed$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FRed$Q10, stand_FRed$Q10, stand_FRed$Q90, stand_FRed$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(4/4*(b-a)+a), legend=("Func. redundancy"), bty="n", cex=1.3)
legend(x=5, y=(4/4*(b-a)+a), legend=expression(paste("+0.19% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Functional Richness #####
FRic <- subset(response_stan, Response == "FRic")
FRic <- FRic$estimate[!is.na(FRic$estimate)]
#average functional richness= 1372.61497
percChange_perYr<-(FRic*100)
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n')
title(xlab="Annual % change", line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="FRic", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
FRic_Est <- subset(Ests, Response == "func_rich")
stand_FRic <- lapply(FRic_Est[,2:11],"*", 100)
yy <- (2.2/4*(b-a)+a)
points(x=stand_FRic$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FRic$Q2.5, stand_FRic$Q2.5, stand_FRic$Q97.5, stand_FRic$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FRic$Q5, stand_FRic$Q5, stand_FRic$Q95, stand_FRic$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FRic$Q10, stand_FRic$Q10, stand_FRic$Q90, stand_FRic$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(3/4*(b-a)+a), legend=("Func. richness"), bty="n", cex=1.3)
legend(x=5, y=(3/4*(b-a)+a), legend=expression(paste("+1% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Func evenness #####
FEve <- subset(response_stan, Response == "FEve")
FEve <- FEve$estimate[!is.na(FEve$estimate)]
#average Func even = 0.517910888
ave_FEve <- 0.517910888
percChange_perYr<-(FEve/ave_FEve)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-1
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
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
yy <- (1.2/4*(b-a)+a)
points(x=stand_FEve$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FEve$Q2.5, stand_FEve$Q2.5, stand_FEve$Q97.5, stand_FEve$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FEve$Q5, stand_FEve$Q5, stand_FEve$Q95, stand_FEve$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FEve$Q10, stand_FEve$Q10, stand_FEve$Q90, stand_FEve$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(2/4*(b-a)+a), legend=("Func. evenness"), bty="n", cex=1.3)
legend(x=5, y=(2/4*(b-a)+a), legend=expression(paste("-0.51% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Func turnover #####
F_to <- subset(response_stan, Response == "F_to")
F_to <- F_to$estimate[!is.na(F_to$estimate)]
#average Func turnover = 0.086496081
#ave_F_to <- 0.086496081
percChange_perYr <- ((exp(F_to)-1)*100)
head(percChange_perYr2)
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-15,15),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
#axis(2, at=0, labels="FTurno", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "tomato", border = "tomato", lwd =2)
F_to_Est <- subset(Ests, Response == "func_turnover")
stand_F_to <- (exp(F_to_Est[,2:11])-1)*100
yy <- (0.2/4*(b-a)+a)
points(x=stand_F_to$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_F_to$Q2.5, stand_F_to$Q2.5, stand_F_to$Q97.5, stand_F_to$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_F_to$Q5, stand_F_to$Q5, stand_F_to$Q95, stand_F_to$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_F_to$Q10, stand_F_to$Q10, stand_F_to$Q90, stand_F_to$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-18, y=(1/4*(b-a)+a), legend=("Func. turnover"), bty="n", cex=1.3)
legend(x=5, y=(1/4*(b-a)+a), legend=expression(paste("-1.41% ",y^-1,, sep = "")), bty="n", cex=1.3)

##
box(lwd=2)
abline(v=0, lwd=1, lty=2)
##
########################################