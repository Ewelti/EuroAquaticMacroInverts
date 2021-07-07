##Set working directory
#setwd("C:/Users/ewelti/Desktop/aquatic_data")
#setwd(dir = "C:/Users/Ellen/Desktop/aquatic_data/")

library(scales)

# attach data
go <- read.csv("outputs/All_siteLevel_and_glmOutput.csv", header=T) # change file name according to the time series to be analyzed
attach(go)
head(go)

par(mfrow=c(2,3),mar=c(2,4,0.4,0.4))
goo <- go[!is.na(SppRich_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$SppRich_Est~goo$meanYear)
summary(mol)
plot(goo$SppRich_Est~goo$meanYear,ylab="Slope of species richness", xlab="", type="n", las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$SppRich_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$SppRich_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

goo <- go[!is.na(ShanH_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$ShanH_Est~goo$meanYear)
summary(mol)
plot(goo$ShanH_Est~goo$meanYear,ylab="Slope of Shannon's H", xlab="", type="n", las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$ShanH_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$ShanH_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

goo <- go[!is.na(EvenJ_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$EvenJ_Est~goo$meanYear)
summary(mol)
plot(goo$EvenJ_Est~goo$meanYear,ylab="Slope of evenness (J)",xlab="", type="n", las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$EvenJ_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$EvenJ_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

par(mar=c(4,4,0.4,0.4))
goo <- go[!is.na(TurnO_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$TurnO_Est~goo$meanYear)
summary(mol)
plot(goo$TurnO_Est~goo$meanYear,ylab="Slope of taxonomic temporal turnover",xlab="", type="n", las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$TurnO_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$TurnO_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

goo <- go[!is.na(Adun_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$Adun_Est~goo$meanYear)
summary(mol)
plot(goo$Adun_Est~goo$meanYear,ylab="Slope of Abundance", xlab="Mean sampling year", type="n", las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$Adun_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$Adun_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

#
############################################################################
#Functional 
par(mfrow=c(2,3),mar=c(2,4,0.4,0.4))

goo <- go[!is.na(FRic_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$FRic_Estgoo$~meanYear)
summary(mol)
plot(goo$FRic_Est~goo$meanYear,ylab="Slope of functional richness", xlab="", type="n", las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$FRic_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$FRic_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

go <- go[order(go$meanYear),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
mol<-lm(goo$RaoQ_Est~goo$meanYear)
summary(mol)
plot(goo$RaoQ_Est~goo$meanYear,ylab="Slope of Rao's Q", xlab="", type="n", las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$RaoQ_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$RaoQ_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

goo <- go[!is.na(FEve_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$FEve_Est~goo$meanYear)
summary(mol)
plot(goo$FEve_Est~goo$meanYear,ylab="Slope of functional evenness", xlab="",type="n", las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$FEve_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$FEve_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

par(mar=c(4,4,0.4,0.4))
goo <- go[!is.na(F_to_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$F_to_Est~goo$meanYear)
summary(mol)
plot(goo$F_to_Est~goo$meanYear,ylab="Slope of functional turnover", xlab="",type="n", las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$F_to_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$F_to_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

goo <- go[!is.na(FDiv_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$FDiv_Est~vmeanYear)
summary(mol)
plot(goo$FDiv_Est~goo$meanYear,ylab="Slope of functional divergence", xlab="Mean sampling year", type="n",las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$FDiv_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$FDiv_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

############################################################################3
#Alien 
par(mfrow=c(3,2),mar=c(2,4,0.4,0.4))

goo <- go[!is.na(AlienAbun_Est),]
head(goo)
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$AlienAbun_Est~goo$meanYear)
summary(mol)
plot(goo$AlienAbun_Est~goo$meanYear,ylab="Slope of non-native abundance", xlab="", type="n",las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$AlienAbun_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$AlienAbun_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

goo <- go[!is.na(AlienSppRich_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$AlienSppRich_Est~goo$meanYear)
summary(mol)
plot(goo$AlienSppRich_Est~goo$meanYear,ylab="Slope of non-native species richness", xlab="",type="n",las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$AlienSppRich_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$AlienSppRich_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

goo <- go[!is.na(nativeAbun_Est),]
head(goo)
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$nativeAbun_Est~goo$meanYear)
summary(mol)
plot(goo$nativeAbun_Est~goo$meanYear,ylab="Slope of native abundance", xlab="", type="n",las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$nativeAbun_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$nativeAbun_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

goo <- go[!is.na(nativeSppRich_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$nativeSppRich_Est~goo$meanYear)
summary(mol)
plot(goo$nativeSppRich_Est~goo$meanYear,ylab="Slope of native species richness", xlab="", type="n",las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$nativeSppRich_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$nativeSppRich_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

par(mar=c(4,4,0.4,0.4))
goo <- go[!is.na(AlienAbunProp_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$AlienAbunProp_Est~goo$meanYear)
summary(mol)
plot(goo$AlienAbunProp_Est~goo$meanYear,ylab="Slope of proportion non-native:native abundance", xlab="Mean sampling year",type="n",las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$AlienAbunProp_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$AlienAbunProp_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)

goo <- go[!is.na(AlienSppRichProp_Est),]
goo$meanYear<-(goo$Ending_year+goo$Starting_year)/2
goo <- goo[order(goo$meanYear),]
mol<-lm(goo$AlienSppRichProp_Est~goo$meanYear)
summary(mol)
plot(goo$AlienSppRichProp_Est~goo$meanYear,ylab="Slope of proportion non-native:native species richness", xlab="Mean sampling year",type="n",las=1)
polygon(x = c(0, 0, 2040, 2040), y = c(-10, 0, 0, -10), col = "grey80",border = "grey80") #abline(h=0,lty=2)
points(goo$AlienSppRichProp_Est~goo$meanYear,pch=19, col=alpha(1,0.2))
box()
plx<-predict(loess(goo$AlienSppRichProp_Est~goo$meanYear), se=T)
lines(goo$meanYear,plx$fit, col="red", lwd=3)
lines(goo$meanYear,plx$fit - qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
lines(goo$meanYear,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, col="red", lwd=1)
##
