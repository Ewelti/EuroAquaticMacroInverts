##Set working directory
setwd("C:/Users/ewelti/Desktop/aquatic_data/movingWindow")
#setwd(dir = "C:/Users/Ellen/Desktop/aquatic_data/movingWindow")


# attach data
go <- read.csv("NorthSouthMovingWindow.csv", header=T) # change file name according to the time series to be analyzed
head(go)

par(mfrow=c(2,3),mar=c(2,5,0.4,0.4))
##
gg <- go[which(go$site_count > 49),]
goo <- gg[!is.na(gg$sppRich),]
head(goo)

nor <- goo[ which(goo$Lat50N=='north'), ]
head(nor)
nor <- nor[order(nor$meanyear),]

sou <- goo[ which(goo$Lat50N=='south'), ]
head(sou)
sou <- sou[order(sou$meanyear),]
######################################################################################
##sppRich
plot(goo$sppRich~goo$meanyear,ylab="", xlab="",
cex=1.5, type="n", las=1, ylim=c(-1.3,1))
title(xlab="Mean year of moving window", line=3,cex.lab=1.5)
title(ylab="Slope of species richness +/- SE", line=3,cex.lab=1.5)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
box()
#
nor$se_sppRich<-as.numeric(nor$se_sppRich)
points(nor$sppRich~nor$meanyear,pch=17, col="darkblue",cex=1.3)
points(nor$sppRich~nor$meanyear,type="l",lwd=2, col="darkblue")
arrows(nor$meanyear, nor$sppRich-nor$se_sppRich, 
nor$meanyear, nor$sppRich+nor$se_sppRich, length=0.05, angle=90, code=3,lwd=1,col="darkblue")
####
sou$se_sppRich<-as.numeric(sou$se_sppRich)
points(sou$sppRich~sou$meanyear,pch=16, col="tomato",cex=1.3)
points(sou$sppRich~sou$meanyear,type="l",lwd=2,col="tomato")
arrows(sou$meanyear, sou$sppRich-sou$se_sppRich, 
sou$meanyear, sou$sppRich+sou$se_sppRich, length=0.05, angle=90, code=3,lwd=1.5,col="tomato")
#
legend("bottomleft", legend=c("North","South"),pch=c(17,16),col=c("darkblue","tomato"),bty="n",cex=2)

######################################################################################
##abun
plot(goo$abun~goo$meanyear,ylab="", xlab="",
cex=1.5, type="n", las=1, ylim=c(-0.05,0.05))
title(xlab="Mean year of moving window", line=3,cex.lab=1.5)
title(ylab="Slope of abundance +/- SE", line=3.3,cex.lab=1.5)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
box()
#
nor$se_abun<-as.numeric(nor$se_abun)
points(nor$abun~nor$meanyear,pch=17, col="darkblue",cex=1.3)
points(nor$abun~nor$meanyear,type="l",lwd=2, col="darkblue")
arrows(nor$meanyear, nor$abun-nor$se_abun, 
nor$meanyear, nor$abun+nor$se_abun, length=0.05, angle=90, code=3,lwd=1,col="darkblue")
####
sou$se_abun<-as.numeric(sou$se_abun)
points(sou$abun~sou$meanyear,pch=16, col="tomato",cex=1.3)
points(sou$abun~sou$meanyear,type="l",lwd=2,col="tomato")
arrows(sou$meanyear, sou$abun-sou$se_abun, 
sou$meanyear, sou$abun+sou$se_abun, length=0.05, angle=90, code=3,lwd=1.5,col="tomato")
#
#legend("bottomleft", legend=c("North","South"),pch=c(17,16),col=c("darkblue","tomato"),bty="n",cex=2)
######################################################################################
##turnover
plot(goo$turnover~goo$meanyear,ylab="", xlab="",
cex=1.5, type="n", las=1, ylim=c(-0.06,0.02))
title(xlab="Mean year of moving window", line=2.7,cex.lab=1.5)
title(ylab="Slope of turnover +/- SE", line=3.3,cex.lab=1.5)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
box()
#
nor$se_turnover<-as.numeric(nor$se_turnover)
points(nor$turnover~nor$meanyear,pch=17, col="darkblue",cex=1.3)
points(nor$turnover~nor$meanyear,type="l",lwd=2, col="darkblue")
arrows(nor$meanyear, nor$turnover-nor$se_turnover, 
nor$meanyear, nor$turnover+nor$se_turnover, length=0.05, angle=90, code=3,lwd=1,col="darkblue")
####
sou$se_turnover<-as.numeric(sou$se_turnover)
points(sou$turnover~sou$meanyear,pch=16, col="tomato",cex=1.3)
points(sou$turnover~sou$meanyear,type="l",lwd=2,col="tomato")
arrows(sou$meanyear, sou$turnover-sou$se_turnover, 
sou$meanyear, sou$turnover+sou$se_turnover, length=0.05, angle=90, code=3,lwd=1.5,col="tomato")
#
#legend("bottomleft", legend=c("North","South"),pch=c(17,16),col=c("darkblue","tomato"),bty="n",cex=2)
#####################################################################################
par(mar=c(4,5,0.4,0.4))
##sppRichRare
plot(goo$sppRichRare~goo$meanyear,ylab="", xlab="",
cex=1.5, type="n", las=1, ylim=c(-0.85,0.5))
title(xlab="Mean year of moving window", line=2.7,cex.lab=1.5)
title(ylab="Slope of rarefied species richness +/- SE", line=3,cex.lab=1.5)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
box()
#
nor$se_sppRichRare<-as.numeric(nor$se_sppRichRare)
points(nor$sppRichRare~nor$meanyear,pch=17, col="darkblue",cex=1.3)
points(nor$sppRichRare~nor$meanyear,type="l",lwd=2, col="darkblue")
arrows(nor$meanyear, nor$sppRichRare-nor$se_sppRichRare, 
nor$meanyear, nor$sppRichRare+nor$se_sppRichRare, length=0.05, angle=90, code=3,lwd=1,col="darkblue")
####
sou$se_sppRichRare<-as.numeric(sou$se_sppRichRare)
points(sou$sppRichRare~sou$meanyear,pch=16, col="tomato",cex=1.3)
points(sou$sppRichRare~sou$meanyear,type="l",lwd=2,col="tomato")
arrows(sou$meanyear, sou$sppRichRare-sou$se_sppRichRare, 
sou$meanyear, sou$sppRichRare+sou$se_sppRichRare, length=0.05, angle=90, code=3,lwd=1.5,col="tomato")
#
#legend("bottomleft", legend=c("North","South"),pch=c(17,16),col=c("darkblue","tomato"),bty="n",cex=2)
################################################################################
##S_PIE
plot(goo$S_PIE~goo$meanyear,ylab="", xlab="",
cex=1.5, type="n", las=1, ylim=c(-0.015,0.028))
title(xlab="Mean year of moving window", line=2.7,cex.lab=1.5)
title(ylab="Slope of S PIE +/- SE", line=3.3,cex.lab=1.5)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
box()
#
nor$se_S_PIE<-as.numeric(nor$se_S_PIE)
points(nor$S_PIE~nor$meanyear,pch=17, col="darkblue",cex=1.3)
points(nor$S_PIE~nor$meanyear,type="l",lwd=2, col="darkblue")
arrows(nor$meanyear, nor$S_PIE-nor$se_S_PIE, 
nor$meanyear, nor$S_PIE+nor$se_S_PIE, length=0.05, angle=90, code=3,lwd=1,col="darkblue")
####
sou$se_S_PIE<-as.numeric(sou$se_S_PIE)
points(sou$S_PIE~sou$meanyear,pch=16, col="tomato",cex=1.3)
points(sou$S_PIE~sou$meanyear,type="l",lwd=2,col="tomato")
arrows(sou$meanyear, sou$S_PIE-sou$se_S_PIE, 
sou$meanyear, sou$S_PIE+sou$se_S_PIE, length=0.05, angle=90, code=3,lwd=1.5,col="tomato")
#
#legend("bottomleft", legend=c("North","South"),pch=c(17,16),col=c("darkblue","tomato"),bty="n",cex=2)
#####################################################################################
##ShanH
plot(goo$ShanH~goo$meanyear,ylab="", xlab="",
cex=1.5, type="n", las=1, ylim=c(-0.03,0.06))
title(xlab="Mean year of moving window", line=2.7,cex.lab=1.5)
title(ylab="Slope of Shannon's diversity +/- SE", line=3.3,cex.lab=1.5)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
box()
#
nor$se_ShanH<-as.numeric(nor$se_ShanH)
points(nor$ShanH~nor$meanyear,pch=17, col="darkblue",cex=1.3)
points(nor$ShanH~nor$meanyear,type="l",lwd=2, col="darkblue")
arrows(nor$meanyear, nor$ShanH-nor$se_ShanH, 
nor$meanyear, nor$ShanH+nor$se_ShanH, length=0.05, angle=90, code=3,lwd=1,col="darkblue")
####
sou$se_ShanH<-as.numeric(sou$se_ShanH)
points(sou$ShanH~sou$meanyear,pch=16, col="tomato",cex=1.3)
points(sou$ShanH~sou$meanyear,type="l",lwd=2,col="tomato")
arrows(sou$meanyear, sou$ShanH-sou$se_ShanH, 
sou$meanyear, sou$ShanH+sou$se_ShanH, length=0.05, angle=90, code=3,lwd=1.5,col="tomato")
#
#legend("bottomleft", legend=c("North","South"),pch=c(17,16),col=c("darkblue","tomato"),bty="n",cex=2)
##