##############################################
#### driver marginal effects ####
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")
meta <- read.csv("outputs/All_siteLevel_and_glmOutput.csv")
driv_m <- median(meta$crop_meanPerc_upstr, na.rm = TRUE)
driv_sd <- sd(meta$crop_meanPerc_upstr, na.rm = TRUE)

tiff(filename = "plots/MarginalEffects/Marg_main_crop.tiff", width = 10, height = 5, units = 'in', res = 600, compression = 'lzw')

##layout
layout(mat = matrix(c(1:15), 
                        nrow = 3, 
                        ncol = 5),
       heights = c(2,2,0.4,2,2,0.4,2,2,0.4,2,2,0.4,2,2,0.4),
       widths = c(0.6,2,2,2,2,0.6,2,2,2,2,0.6,2,2,2,2))

par(mar=c(2,3.8,0.4,0.8))
#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab="Predicted trend", line=-1.2, cex.lab=1.6)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab="Predicted trend", line=-1.2, cex.lab=1.6)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

##richness
ma <- readRDS("outputs/Driver_marginalEffects/preds-scrop_meanPerc_upstr-spp_richness.rds")
head(ma)
driv <- (ma$scrop_meanPerc_upstr*driv_sd+driv_m)*100
plot(1, type="n", xlim=c(min(driv),max(driv)), ylim=c(min(ma$lower),max(ma$upper)),xaxs="i",las=1,ylab="",xlab="")
polygon(x = c(-1000, -1000, 1000, 1000), y = c(-1000, 0, 0, -1000), col = "mistyrose1", lty = 0)   
polygon(c(driv, rev(driv)), c(ma$lower, rev(ma$upper)), col = "lightgrey", lty = 0)
abline(h=0, lty=2,col="red")
points(ma$trend ~ driv,pch=21,type="l",cex=2,lwd=4)
legend("topleft", legend=("a, Taxon richness"), bty="n", cex=1.3)
box(lwd=2)

##func richness
ma <- readRDS("outputs/Driver_marginalEffects/preds-scrop_meanPerc_upstr-FRic.rds")
driv <- (ma$scrop_meanPerc_upstr*driv_sd+driv_m)*100
plot(1, type="n", xlim=c(min(driv),max(driv)), ylim=c(min(ma$lower),max(ma$upper)),xaxs="i",las=1,ylab="",xlab="")
polygon(x = c(-1000, -1000, 1000, 1000), y = c(-1000, 0, 0, -1000), col = "mistyrose1", lty = 0)   
polygon(c(driv, rev(driv)), c(ma$lower, rev(ma$upper)), col = "lightgrey", lty = 0)
abline(h=0, lty=2,col="red")
points(ma$trend ~ driv,pch=21,type="l",cex=2,lwd=4)
legend("topleft", legend=("e, Func. richness"), bty="n", cex=1.3)
box(lwd=2)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="% area crop upstream", line=-0.5, cex.lab=1.4)

##abundance
ma <- readRDS("outputs/Driver_marginalEffects/preds-scrop_meanPerc_upstr-abundance.rds")
driv <- (ma$scrop_meanPerc_upstr*driv_sd+driv_m)*100
plot(1, type="n", xlim=c(min(driv),max(driv)), ylim=c(min(ma$lower),max(ma$upper)),xaxs="i",las=1,ylab="",xlab="")
polygon(x = c(-1000, -1000, 1000, 1000), y = c(-1000, 0, 0, -1000), col = "mistyrose1", lty = 0)   
polygon(c(driv, rev(driv)), c(ma$lower, rev(ma$upper)), col = "lightgrey", lty = 0)
abline(h=0, lty=2,col="red")
points(ma$trend ~ driv,pch=21,type="l",cex=2,lwd=4)
legend("topleft", legend=("b, Abundance"), bty="n", cex=1.3)
box(lwd=2)

##functional redundancy
ma <- readRDS("outputs/Driver_marginalEffects/preds-scrop_meanPerc_upstr-FRed.rds")
driv <- (ma$scrop_meanPerc_upstr*driv_sd+driv_m)*100
plot(1, type="n", xlim=c(min(driv),max(driv)), ylim=c(min(ma$lower),max(ma$upper)),xaxs="i",las=1,ylab="",xlab="")
polygon(x = c(-1000, -1000, 1000, 1000), y = c(-1000, 0, 0, -1000), col = "mistyrose1", lty = 0)   
polygon(c(driv, rev(driv)), c(ma$lower, rev(ma$upper)), col = "lightgrey", lty = 0)
abline(h=0, lty=2,col="red")
points(ma$trend ~ driv,pch=21,type="l",cex=2,lwd=4)
legend("topleft", legend=("f, Func. redundancy"), bty="n", cex=1.3)
box(lwd=2)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="% area crop upstream", line=-0.5, cex.lab=1.4)

##evenness
ma <- readRDS("outputs/Driver_marginalEffects/preds-scrop_meanPerc_upstr-E10.rds")
driv <- (ma$scrop_meanPerc_upstr*driv_sd+driv_m)*100
plot(1, type="n", xlim=c(min(driv),max(driv)), ylim=c(min(ma$lower),max(ma$upper)),xaxs="i",las=1,ylab="",xlab="")
polygon(x = c(-1000, -1000, 1000, 1000), y = c(-1000, 0, 0, -1000), col = "mistyrose1", lty = 0)   
polygon(c(driv, rev(driv)), c(ma$lower, rev(ma$upper)), col = "lightgrey", lty = 0)
abline(h=0, lty=2,col="red")
points(ma$trend ~ driv,pch=21,type="l",cex=2,lwd=4)
legend("topleft", legend=("c, Evenness"), bty="n", cex=1.3)
box(lwd=2)

##func evenness
ma <- readRDS("outputs/Driver_marginalEffects/preds-scrop_meanPerc_upstr-FEve.rds")
driv <- (ma$scrop_meanPerc_upstr*driv_sd+driv_m)*100
plot(1, type="n", xlim=c(min(driv),max(driv)), ylim=c(min(ma$lower),max(ma$upper)),xaxs="i",las=1,ylab="",xlab="")
polygon(x = c(-1000, -1000, 1000, 1000), y = c(-1000, 0, 0, -1000), col = "mistyrose1", lty = 0)   
polygon(c(driv, rev(driv)), c(ma$lower, rev(ma$upper)), col = "lightgrey", lty = 0)
abline(h=0, lty=2,col="red")
points(ma$trend ~ driv,pch=21,type="l",cex=2,lwd=4)
legend("topleft", legend=("g, Func. evenness"), bty="n", cex=1.3)
box(lwd=2)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="% area crop upstream", line=-0.5, cex.lab=1.4)

##turnover
ma <- readRDS("outputs/Driver_marginalEffects/preds-scrop_meanPerc_upstr-turnover.rds")
driv <- (ma$scrop_meanPerc_upstr*driv_sd+driv_m)*100
plot(1, type="n", xlim=c(min(driv),max(driv)), ylim=c(min(ma$lower),max(ma$upper)),xaxs="i",las=1,ylab="",xlab="")
polygon(x = c(-1000, -1000, 1000, 1000), y = c(-1000, 0, 0, -1000), col = "mistyrose1", lty = 0)   
polygon(c(driv, rev(driv)), c(ma$lower, rev(ma$upper)), col = "lightgrey", lty = 0)
abline(h=0, lty=2,col="red")
points(ma$trend ~ driv,pch=21,type="l",cex=2,lwd=4)
legend("topleft", legend=("d, Turnover"), bty="n", cex=1.3)
box(lwd=2)

##func turnover
ma <- readRDS("outputs/Driver_marginalEffects/preds-scrop_meanPerc_upstr-F_to.rds")
driv <- (ma$scrop_meanPerc_upstr*driv_sd+driv_m)*100
plot(1, type="n", xlim=c(min(driv),max(driv)), ylim=c(min(ma$lower),max(ma$upper)),xaxs="i",las=1,ylab="",xlab="")
polygon(x = c(-1000, -1000, 1000, 1000), y = c(-1000, 0, 0, -1000), col = "mistyrose1", lty = 0)   
polygon(c(driv, rev(driv)), c(ma$lower, rev(ma$upper)), col = "lightgrey", lty = 0)
abline(h=0, lty=2,col="red")
points(ma$trend ~ driv,pch=21,type="l",cex=2,lwd=4)
legend("topleft", legend=("h, Func. turnover"), bty="n", cex=1.3)
box(lwd=2)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="% area crop upstream", line=-0.5, cex.lab=1.4)

##
dev.off()
##

