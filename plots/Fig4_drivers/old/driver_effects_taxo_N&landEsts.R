##Set working directory
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/")

# attach data
#drivEst <- read.csv("outputs/Driver_metaanaly_Ests.csv")
drivEst <- read.csv("outputs/Driver_horseshoe_metaanaly_Ests.csv")
head(drivEst)
unique(drivEst$Response)

options(scipen=999)

#tiff(filename = "plots/drivers/Drivers_taxonDiv.tiff", width = 12, height = 9, units = 'in', res = 600, compression = 'lzw')
tiff(filename = "plots/drivers/Drivers_horseshoe_taxonDiv.tiff", width = 12, height = 9, units = 'in', res = 600, compression = 'lzw')


layout(mat = matrix(c(1:8), 
                        nrow = 2, 
                        ncol = 4),
       heights = c(2,2,2,2,2,2,2,2),
       widths = c(1,2,2,2,1,2,2,2))
#layout.show(8)

par(mar=c(4,0.4,0.4,0.4))
#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
#####################Taxonomic metrics ##############################
#### Spp Richness #####
sub <- subset(drivEst, Response == "spp_richness")
est <- sub$Estimate[2:16]
yy <- c(1,3,2,4,5,6,7,8,9,10,14,12,11,13,15)
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.10,0.12), cex=2)
axis(2, at=yy, labels=c("ppt sl.", "tmax sl.", "ppt mean", "tmax mean", "str. order",
    "accum.", "elevation", "slope", "N sl.", "N mean", "urban mean", "crop mean", "crop sl.",
    "urban sl.", "dam"), las=1,cex.axis=1.3)
title(xlab="Est. on species richness", line=2.4,cex.lab=1.3)
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
box(lwd=2)
points(yy ~ est, pch="l",cex=2)
yyy1=c(0.975,1.025,1.025,0.975)
yyy2=c(0.95,1.05,1.05,0.95)
yyy3=c(0.9,1.1,1.1,0.9)
polygon(x=c(sub$Q2.5[2], sub$Q2.5[2], sub$Q97.5[2], sub$Q97.5[2]),y=yyy1, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[2], sub$Q5[2], sub$Q95[2], sub$Q95[2]),y=yyy2, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[2], sub$Q10[2], sub$Q90[2], sub$Q90[2]),y=yyy3, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[4], sub$Q2.5[4], sub$Q97.5[4], sub$Q97.5[4]),y=(yyy1+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[4], sub$Q5[4], sub$Q95[4], sub$Q95[4]),y=(yyy2+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[4], sub$Q10[4], sub$Q90[4], sub$Q90[4]),y=(yyy3+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[3], sub$Q2.5[3], sub$Q97.5[3], sub$Q97.5[3]),y=(yyy1+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[3], sub$Q5[3], sub$Q95[3], sub$Q95[3]),y=(yyy2+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[3], sub$Q10[3], sub$Q90[3], sub$Q90[3]),y=(yyy3+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[5], sub$Q2.5[5], sub$Q97.5[5], sub$Q97.5[5]),y=(yyy1+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[5], sub$Q5[5], sub$Q95[5], sub$Q95[5]),y=(yyy2+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[5], sub$Q10[5], sub$Q90[5], sub$Q90[5]),y=(yyy3+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[6], sub$Q2.5[6], sub$Q97.5[6], sub$Q97.5[6]),y=(yyy1+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[6], sub$Q5[6], sub$Q95[6], sub$Q95[6]),y=(yyy2+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[6], sub$Q10[6], sub$Q90[6], sub$Q90[6]),y=(yyy3+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[7], sub$Q2.5[7], sub$Q97.5[7], sub$Q97.5[7]),y=(yyy1+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[7], sub$Q5[7], sub$Q95[7], sub$Q95[7]),y=(yyy2+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[7], sub$Q10[7], sub$Q90[7], sub$Q90[7]),y=(yyy3+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[8], sub$Q2.5[8], sub$Q97.5[8], sub$Q97.5[8]),y=(yyy1+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[8], sub$Q5[8], sub$Q95[8], sub$Q95[8]),y=(yyy2+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[8], sub$Q10[8], sub$Q90[8], sub$Q90[8]),y=(yyy3+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[9], sub$Q2.5[9], sub$Q97.5[9], sub$Q97.5[9]),y=(yyy1+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[9], sub$Q5[9], sub$Q95[9], sub$Q95[9]),y=(yyy2+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[9], sub$Q10[9], sub$Q90[9], sub$Q90[9]),y=(yyy3+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[10], sub$Q2.5[10], sub$Q97.5[10], sub$Q97.5[10]),y=(yyy1+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[10], sub$Q5[10], sub$Q95[10], sub$Q95[10]),y=(yyy2+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[10], sub$Q10[10], sub$Q90[10], sub$Q90[10]),y=(yyy3+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[11], sub$Q2.5[11], sub$Q97.5[11], sub$Q97.5[11]),y=(yyy1+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[11], sub$Q5[11], sub$Q95[11], sub$Q95[11]),y=(yyy2+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[11], sub$Q10[11], sub$Q90[11], sub$Q90[11]),y=(yyy3+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[14], sub$Q2.5[14], sub$Q97.5[14], sub$Q97.5[14]),y=(yyy1+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[14], sub$Q5[14], sub$Q95[14], sub$Q95[14]),y=(yyy2+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[14], sub$Q10[14], sub$Q90[14], sub$Q90[14]),y=(yyy3+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[13], sub$Q2.5[13], sub$Q97.5[13], sub$Q97.5[13]),y=(yyy1+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[13], sub$Q5[13], sub$Q95[13], sub$Q95[13]),y=(yyy2+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[13], sub$Q10[13], sub$Q90[13], sub$Q90[13]),y=(yyy3+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[15], sub$Q2.5[15], sub$Q97.5[15], sub$Q97.5[15]),y=(yyy1+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[15], sub$Q5[15], sub$Q95[15], sub$Q95[15]),y=(yyy2+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[15], sub$Q10[15], sub$Q90[15], sub$Q90[15]),y=(yyy3+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[12], sub$Q2.5[12], sub$Q97.5[12], sub$Q97.5[12]),y=(yyy1+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[12], sub$Q5[12], sub$Q95[12], sub$Q95[12]),y=(yyy2+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[12], sub$Q10[12], sub$Q90[12], sub$Q90[12]),y=(yyy3+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[16], sub$Q2.5[16], sub$Q97.5[16], sub$Q97.5[16]),y=(yyy1+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[16], sub$Q5[16], sub$Q95[16], sub$Q95[16]),y=(yyy2+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[16], sub$Q10[16], sub$Q90[16], sub$Q90[16]),y=(yyy3+14), col = 1, border = 1,lwd = 1)
abline(h=4.5,lty=2)
abline(h=8.5,lty=2)
abline(h=10.5,lty=2)
abline(h=14.5,lty=2)

#### shannonsH #####
sub <- subset(drivEst, Response == "shannonsH")
est <- sub$Estimate[2:16]
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.005,0.007), cex=2)
title(xlab="Est. on Shannon's H", line=2.4,cex.lab=1.3)
axis(2, at=yy, labels=c("ppt sl.", "tmax sl.", "ppt mean", "tmax mean", "str. order",
    "accum.", "elevation", "slope", "N sl.", "N mean", "urban mean", "crop mean", "crop sl.",
    "urban sl.", "dam"), las=1,cex.axis=1.3)
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
box(lwd=2)
points(yy ~ est, pch="l",cex=2)
polygon(x=c(sub$Q2.5[2], sub$Q2.5[2], sub$Q97.5[2], sub$Q97.5[2]),y=yyy1, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[2], sub$Q5[2], sub$Q95[2], sub$Q95[2]),y=yyy2, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[2], sub$Q10[2], sub$Q90[2], sub$Q90[2]),y=yyy3, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[4], sub$Q2.5[4], sub$Q97.5[4], sub$Q97.5[4]),y=(yyy1+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[4], sub$Q5[4], sub$Q95[4], sub$Q95[4]),y=(yyy2+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[4], sub$Q10[4], sub$Q90[4], sub$Q90[4]),y=(yyy3+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[3], sub$Q2.5[3], sub$Q97.5[3], sub$Q97.5[3]),y=(yyy1+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[3], sub$Q5[3], sub$Q95[3], sub$Q95[3]),y=(yyy2+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[3], sub$Q10[3], sub$Q90[3], sub$Q90[3]),y=(yyy3+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[5], sub$Q2.5[5], sub$Q97.5[5], sub$Q97.5[5]),y=(yyy1+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[5], sub$Q5[5], sub$Q95[5], sub$Q95[5]),y=(yyy2+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[5], sub$Q10[5], sub$Q90[5], sub$Q90[5]),y=(yyy3+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[6], sub$Q2.5[6], sub$Q97.5[6], sub$Q97.5[6]),y=(yyy1+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[6], sub$Q5[6], sub$Q95[6], sub$Q95[6]),y=(yyy2+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[6], sub$Q10[6], sub$Q90[6], sub$Q90[6]),y=(yyy3+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[7], sub$Q2.5[7], sub$Q97.5[7], sub$Q97.5[7]),y=(yyy1+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[7], sub$Q5[7], sub$Q95[7], sub$Q95[7]),y=(yyy2+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[7], sub$Q10[7], sub$Q90[7], sub$Q90[7]),y=(yyy3+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[8], sub$Q2.5[8], sub$Q97.5[8], sub$Q97.5[8]),y=(yyy1+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[8], sub$Q5[8], sub$Q95[8], sub$Q95[8]),y=(yyy2+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[8], sub$Q10[8], sub$Q90[8], sub$Q90[8]),y=(yyy3+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[9], sub$Q2.5[9], sub$Q97.5[9], sub$Q97.5[9]),y=(yyy1+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[9], sub$Q5[9], sub$Q95[9], sub$Q95[9]),y=(yyy2+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[9], sub$Q10[9], sub$Q90[9], sub$Q90[9]),y=(yyy3+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[10], sub$Q2.5[10], sub$Q97.5[10], sub$Q97.5[10]),y=(yyy1+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[10], sub$Q5[10], sub$Q95[10], sub$Q95[10]),y=(yyy2+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[10], sub$Q10[10], sub$Q90[10], sub$Q90[10]),y=(yyy3+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[11], sub$Q2.5[11], sub$Q97.5[11], sub$Q97.5[11]),y=(yyy1+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[11], sub$Q5[11], sub$Q95[11], sub$Q95[11]),y=(yyy2+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[11], sub$Q10[11], sub$Q90[11], sub$Q90[11]),y=(yyy3+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[14], sub$Q2.5[14], sub$Q97.5[14], sub$Q97.5[14]),y=(yyy1+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[14], sub$Q5[14], sub$Q95[14], sub$Q95[14]),y=(yyy2+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[14], sub$Q10[14], sub$Q90[14], sub$Q90[14]),y=(yyy3+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[13], sub$Q2.5[13], sub$Q97.5[13], sub$Q97.5[13]),y=(yyy1+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[13], sub$Q5[13], sub$Q95[13], sub$Q95[13]),y=(yyy2+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[13], sub$Q10[13], sub$Q90[13], sub$Q90[13]),y=(yyy3+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[15], sub$Q2.5[15], sub$Q97.5[15], sub$Q97.5[15]),y=(yyy1+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[15], sub$Q5[15], sub$Q95[15], sub$Q95[15]),y=(yyy2+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[15], sub$Q10[15], sub$Q90[15], sub$Q90[15]),y=(yyy3+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[12], sub$Q2.5[12], sub$Q97.5[12], sub$Q97.5[12]),y=(yyy1+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[12], sub$Q5[12], sub$Q95[12], sub$Q95[12]),y=(yyy2+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[12], sub$Q10[12], sub$Q90[12], sub$Q90[12]),y=(yyy3+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[16], sub$Q2.5[16], sub$Q97.5[16], sub$Q97.5[16]),y=(yyy1+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[16], sub$Q5[16], sub$Q95[16], sub$Q95[16]),y=(yyy2+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[16], sub$Q10[16], sub$Q90[16], sub$Q90[16]),y=(yyy3+14), col = 1, border = 1,lwd = 1)
abline(h=4.5,lty=2)
abline(h=8.5,lty=2)
abline(h=10.5,lty=2)
abline(h=14.5,lty=2)

#### spp_richness_rarefied #####
sub <- subset(drivEst, Response == "spp_richness_rarefied")
est <- sub$Estimate[2:16]
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.07,0.05), cex=2)
title(xlab="Est. on rarefied spp richness", line=2.4,cex.lab=1.3)
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
box(lwd=2)
points(yy ~ est, pch="l",cex=2)
polygon(x=c(sub$Q2.5[2], sub$Q2.5[2], sub$Q97.5[2], sub$Q97.5[2]),y=yyy1, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[2], sub$Q5[2], sub$Q95[2], sub$Q95[2]),y=yyy2, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[2], sub$Q10[2], sub$Q90[2], sub$Q90[2]),y=yyy3, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[4], sub$Q2.5[4], sub$Q97.5[4], sub$Q97.5[4]),y=(yyy1+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[4], sub$Q5[4], sub$Q95[4], sub$Q95[4]),y=(yyy2+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[4], sub$Q10[4], sub$Q90[4], sub$Q90[4]),y=(yyy3+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[3], sub$Q2.5[3], sub$Q97.5[3], sub$Q97.5[3]),y=(yyy1+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[3], sub$Q5[3], sub$Q95[3], sub$Q95[3]),y=(yyy2+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[3], sub$Q10[3], sub$Q90[3], sub$Q90[3]),y=(yyy3+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[5], sub$Q2.5[5], sub$Q97.5[5], sub$Q97.5[5]),y=(yyy1+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[5], sub$Q5[5], sub$Q95[5], sub$Q95[5]),y=(yyy2+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[5], sub$Q10[5], sub$Q90[5], sub$Q90[5]),y=(yyy3+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[6], sub$Q2.5[6], sub$Q97.5[6], sub$Q97.5[6]),y=(yyy1+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[6], sub$Q5[6], sub$Q95[6], sub$Q95[6]),y=(yyy2+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[6], sub$Q10[6], sub$Q90[6], sub$Q90[6]),y=(yyy3+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[7], sub$Q2.5[7], sub$Q97.5[7], sub$Q97.5[7]),y=(yyy1+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[7], sub$Q5[7], sub$Q95[7], sub$Q95[7]),y=(yyy2+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[7], sub$Q10[7], sub$Q90[7], sub$Q90[7]),y=(yyy3+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[8], sub$Q2.5[8], sub$Q97.5[8], sub$Q97.5[8]),y=(yyy1+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[8], sub$Q5[8], sub$Q95[8], sub$Q95[8]),y=(yyy2+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[8], sub$Q10[8], sub$Q90[8], sub$Q90[8]),y=(yyy3+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[9], sub$Q2.5[9], sub$Q97.5[9], sub$Q97.5[9]),y=(yyy1+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[9], sub$Q5[9], sub$Q95[9], sub$Q95[9]),y=(yyy2+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[9], sub$Q10[9], sub$Q90[9], sub$Q90[9]),y=(yyy3+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[10], sub$Q2.5[10], sub$Q97.5[10], sub$Q97.5[10]),y=(yyy1+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[10], sub$Q5[10], sub$Q95[10], sub$Q95[10]),y=(yyy2+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[10], sub$Q10[10], sub$Q90[10], sub$Q90[10]),y=(yyy3+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[11], sub$Q2.5[11], sub$Q97.5[11], sub$Q97.5[11]),y=(yyy1+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[11], sub$Q5[11], sub$Q95[11], sub$Q95[11]),y=(yyy2+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[11], sub$Q10[11], sub$Q90[11], sub$Q90[11]),y=(yyy3+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[14], sub$Q2.5[14], sub$Q97.5[14], sub$Q97.5[14]),y=(yyy1+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[14], sub$Q5[14], sub$Q95[14], sub$Q95[14]),y=(yyy2+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[14], sub$Q10[14], sub$Q90[14], sub$Q90[14]),y=(yyy3+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[13], sub$Q2.5[13], sub$Q97.5[13], sub$Q97.5[13]),y=(yyy1+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[13], sub$Q5[13], sub$Q95[13], sub$Q95[13]),y=(yyy2+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[13], sub$Q10[13], sub$Q90[13], sub$Q90[13]),y=(yyy3+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[15], sub$Q2.5[15], sub$Q97.5[15], sub$Q97.5[15]),y=(yyy1+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[15], sub$Q5[15], sub$Q95[15], sub$Q95[15]),y=(yyy2+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[15], sub$Q10[15], sub$Q90[15], sub$Q90[15]),y=(yyy3+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[12], sub$Q2.5[12], sub$Q97.5[12], sub$Q97.5[12]),y=(yyy1+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[12], sub$Q5[12], sub$Q95[12], sub$Q95[12]),y=(yyy2+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[12], sub$Q10[12], sub$Q90[12], sub$Q90[12]),y=(yyy3+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[16], sub$Q2.5[16], sub$Q97.5[16], sub$Q97.5[16]),y=(yyy1+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[16], sub$Q5[16], sub$Q95[16], sub$Q95[16]),y=(yyy2+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[16], sub$Q10[16], sub$Q90[16], sub$Q90[16]),y=(yyy3+14), col = 1, border = 1,lwd = 1)
abline(h=4.5,lty=2)
abline(h=8.5,lty=2)
abline(h=10.5,lty=2)
abline(h=14.5,lty=2)

#### E10 #####
sub <- subset(drivEst, Response == "E10")
est <- sub$Estimate[2:16]
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.0014,0.001), cex=2)
title(xlab="Est. on evenness", line=2.4,cex.lab=1.3)
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
box(lwd=2)
points(yy ~ est, pch="l",cex=2)
polygon(x=c(sub$Q2.5[2], sub$Q2.5[2], sub$Q97.5[2], sub$Q97.5[2]),y=yyy1, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[2], sub$Q5[2], sub$Q95[2], sub$Q95[2]),y=yyy2, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[2], sub$Q10[2], sub$Q90[2], sub$Q90[2]),y=yyy3, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[4], sub$Q2.5[4], sub$Q97.5[4], sub$Q97.5[4]),y=(yyy1+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[4], sub$Q5[4], sub$Q95[4], sub$Q95[4]),y=(yyy2+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[4], sub$Q10[4], sub$Q90[4], sub$Q90[4]),y=(yyy3+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[3], sub$Q2.5[3], sub$Q97.5[3], sub$Q97.5[3]),y=(yyy1+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[3], sub$Q5[3], sub$Q95[3], sub$Q95[3]),y=(yyy2+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[3], sub$Q10[3], sub$Q90[3], sub$Q90[3]),y=(yyy3+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[5], sub$Q2.5[5], sub$Q97.5[5], sub$Q97.5[5]),y=(yyy1+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[5], sub$Q5[5], sub$Q95[5], sub$Q95[5]),y=(yyy2+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[5], sub$Q10[5], sub$Q90[5], sub$Q90[5]),y=(yyy3+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[6], sub$Q2.5[6], sub$Q97.5[6], sub$Q97.5[6]),y=(yyy1+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[6], sub$Q5[6], sub$Q95[6], sub$Q95[6]),y=(yyy2+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[6], sub$Q10[6], sub$Q90[6], sub$Q90[6]),y=(yyy3+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[7], sub$Q2.5[7], sub$Q97.5[7], sub$Q97.5[7]),y=(yyy1+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[7], sub$Q5[7], sub$Q95[7], sub$Q95[7]),y=(yyy2+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[7], sub$Q10[7], sub$Q90[7], sub$Q90[7]),y=(yyy3+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[8], sub$Q2.5[8], sub$Q97.5[8], sub$Q97.5[8]),y=(yyy1+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[8], sub$Q5[8], sub$Q95[8], sub$Q95[8]),y=(yyy2+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[8], sub$Q10[8], sub$Q90[8], sub$Q90[8]),y=(yyy3+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[9], sub$Q2.5[9], sub$Q97.5[9], sub$Q97.5[9]),y=(yyy1+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[9], sub$Q5[9], sub$Q95[9], sub$Q95[9]),y=(yyy2+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[9], sub$Q10[9], sub$Q90[9], sub$Q90[9]),y=(yyy3+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[10], sub$Q2.5[10], sub$Q97.5[10], sub$Q97.5[10]),y=(yyy1+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[10], sub$Q5[10], sub$Q95[10], sub$Q95[10]),y=(yyy2+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[10], sub$Q10[10], sub$Q90[10], sub$Q90[10]),y=(yyy3+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[11], sub$Q2.5[11], sub$Q97.5[11], sub$Q97.5[11]),y=(yyy1+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[11], sub$Q5[11], sub$Q95[11], sub$Q95[11]),y=(yyy2+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[11], sub$Q10[11], sub$Q90[11], sub$Q90[11]),y=(yyy3+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[14], sub$Q2.5[14], sub$Q97.5[14], sub$Q97.5[14]),y=(yyy1+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[14], sub$Q5[14], sub$Q95[14], sub$Q95[14]),y=(yyy2+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[14], sub$Q10[14], sub$Q90[14], sub$Q90[14]),y=(yyy3+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[13], sub$Q2.5[13], sub$Q97.5[13], sub$Q97.5[13]),y=(yyy1+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[13], sub$Q5[13], sub$Q95[13], sub$Q95[13]),y=(yyy2+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[13], sub$Q10[13], sub$Q90[13], sub$Q90[13]),y=(yyy3+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[15], sub$Q2.5[15], sub$Q97.5[15], sub$Q97.5[15]),y=(yyy1+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[15], sub$Q5[15], sub$Q95[15], sub$Q95[15]),y=(yyy2+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[15], sub$Q10[15], sub$Q90[15], sub$Q90[15]),y=(yyy3+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[12], sub$Q2.5[12], sub$Q97.5[12], sub$Q97.5[12]),y=(yyy1+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[12], sub$Q5[12], sub$Q95[12], sub$Q95[12]),y=(yyy2+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[12], sub$Q10[12], sub$Q90[12], sub$Q90[12]),y=(yyy3+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[16], sub$Q2.5[16], sub$Q97.5[16], sub$Q97.5[16]),y=(yyy1+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[16], sub$Q5[16], sub$Q95[16], sub$Q95[16]),y=(yyy2+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[16], sub$Q10[16], sub$Q90[16], sub$Q90[16]),y=(yyy3+14), col = 1, border = 1,lwd = 1)
abline(h=4.5,lty=2)
abline(h=8.5,lty=2)
abline(h=10.5,lty=2)
abline(h=14.5,lty=2)

#### Abundance #####
sub <- subset(drivEst, Response == "abundance")
est <- sub$Estimate[2:16]
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.004,0.004), cex=2)
title(xlab="Est. on abundance", line=2.4,cex.lab=1.3)
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
box(lwd=2)
points(yy ~ est, pch="l",cex=2)
polygon(x=c(sub$Q2.5[2], sub$Q2.5[2], sub$Q97.5[2], sub$Q97.5[2]),y=yyy1, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[2], sub$Q5[2], sub$Q95[2], sub$Q95[2]),y=yyy2, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[2], sub$Q10[2], sub$Q90[2], sub$Q90[2]),y=yyy3, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[4], sub$Q2.5[4], sub$Q97.5[4], sub$Q97.5[4]),y=(yyy1+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[4], sub$Q5[4], sub$Q95[4], sub$Q95[4]),y=(yyy2+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[4], sub$Q10[4], sub$Q90[4], sub$Q90[4]),y=(yyy3+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[3], sub$Q2.5[3], sub$Q97.5[3], sub$Q97.5[3]),y=(yyy1+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[3], sub$Q5[3], sub$Q95[3], sub$Q95[3]),y=(yyy2+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[3], sub$Q10[3], sub$Q90[3], sub$Q90[3]),y=(yyy3+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[5], sub$Q2.5[5], sub$Q97.5[5], sub$Q97.5[5]),y=(yyy1+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[5], sub$Q5[5], sub$Q95[5], sub$Q95[5]),y=(yyy2+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[5], sub$Q10[5], sub$Q90[5], sub$Q90[5]),y=(yyy3+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[6], sub$Q2.5[6], sub$Q97.5[6], sub$Q97.5[6]),y=(yyy1+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[6], sub$Q5[6], sub$Q95[6], sub$Q95[6]),y=(yyy2+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[6], sub$Q10[6], sub$Q90[6], sub$Q90[6]),y=(yyy3+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[7], sub$Q2.5[7], sub$Q97.5[7], sub$Q97.5[7]),y=(yyy1+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[7], sub$Q5[7], sub$Q95[7], sub$Q95[7]),y=(yyy2+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[7], sub$Q10[7], sub$Q90[7], sub$Q90[7]),y=(yyy3+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[8], sub$Q2.5[8], sub$Q97.5[8], sub$Q97.5[8]),y=(yyy1+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[8], sub$Q5[8], sub$Q95[8], sub$Q95[8]),y=(yyy2+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[8], sub$Q10[8], sub$Q90[8], sub$Q90[8]),y=(yyy3+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[9], sub$Q2.5[9], sub$Q97.5[9], sub$Q97.5[9]),y=(yyy1+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[9], sub$Q5[9], sub$Q95[9], sub$Q95[9]),y=(yyy2+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[9], sub$Q10[9], sub$Q90[9], sub$Q90[9]),y=(yyy3+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[10], sub$Q2.5[10], sub$Q97.5[10], sub$Q97.5[10]),y=(yyy1+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[10], sub$Q5[10], sub$Q95[10], sub$Q95[10]),y=(yyy2+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[10], sub$Q10[10], sub$Q90[10], sub$Q90[10]),y=(yyy3+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[11], sub$Q2.5[11], sub$Q97.5[11], sub$Q97.5[11]),y=(yyy1+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[11], sub$Q5[11], sub$Q95[11], sub$Q95[11]),y=(yyy2+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[11], sub$Q10[11], sub$Q90[11], sub$Q90[11]),y=(yyy3+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[14], sub$Q2.5[14], sub$Q97.5[14], sub$Q97.5[14]),y=(yyy1+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[14], sub$Q5[14], sub$Q95[14], sub$Q95[14]),y=(yyy2+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[14], sub$Q10[14], sub$Q90[14], sub$Q90[14]),y=(yyy3+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[13], sub$Q2.5[13], sub$Q97.5[13], sub$Q97.5[13]),y=(yyy1+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[13], sub$Q5[13], sub$Q95[13], sub$Q95[13]),y=(yyy2+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[13], sub$Q10[13], sub$Q90[13], sub$Q90[13]),y=(yyy3+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[15], sub$Q2.5[15], sub$Q97.5[15], sub$Q97.5[15]),y=(yyy1+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[15], sub$Q5[15], sub$Q95[15], sub$Q95[15]),y=(yyy2+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[15], sub$Q10[15], sub$Q90[15], sub$Q90[15]),y=(yyy3+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[12], sub$Q2.5[12], sub$Q97.5[12], sub$Q97.5[12]),y=(yyy1+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[12], sub$Q5[12], sub$Q95[12], sub$Q95[12]),y=(yyy2+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[12], sub$Q10[12], sub$Q90[12], sub$Q90[12]),y=(yyy3+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[16], sub$Q2.5[16], sub$Q97.5[16], sub$Q97.5[16]),y=(yyy1+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[16], sub$Q5[16], sub$Q95[16], sub$Q95[16]),y=(yyy2+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[16], sub$Q10[16], sub$Q90[16], sub$Q90[16]),y=(yyy3+14), col = 1, border = 1,lwd = 1)
abline(h=4.5,lty=2)
abline(h=8.5,lty=2)
abline(h=10.5,lty=2)
abline(h=14.5,lty=2)

#### Turnover #####
sub <- subset(drivEst, Response == "turnover")
est <- sub$Estimate[2:16]
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.0025,0.0016), cex=2)
title(xlab="Est. on turnover", line=2.4,cex.lab=1.3)
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
box(lwd=2)
points(yy ~ est, pch="l",cex=2)
polygon(x=c(sub$Q2.5[2], sub$Q2.5[2], sub$Q97.5[2], sub$Q97.5[2]),y=yyy1, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[2], sub$Q5[2], sub$Q95[2], sub$Q95[2]),y=yyy2, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[2], sub$Q10[2], sub$Q90[2], sub$Q90[2]),y=yyy3, col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[4], sub$Q2.5[4], sub$Q97.5[4], sub$Q97.5[4]),y=(yyy1+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[4], sub$Q5[4], sub$Q95[4], sub$Q95[4]),y=(yyy2+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[4], sub$Q10[4], sub$Q90[4], sub$Q90[4]),y=(yyy3+1), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[3], sub$Q2.5[3], sub$Q97.5[3], sub$Q97.5[3]),y=(yyy1+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[3], sub$Q5[3], sub$Q95[3], sub$Q95[3]),y=(yyy2+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[3], sub$Q10[3], sub$Q90[3], sub$Q90[3]),y=(yyy3+2), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[5], sub$Q2.5[5], sub$Q97.5[5], sub$Q97.5[5]),y=(yyy1+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[5], sub$Q5[5], sub$Q95[5], sub$Q95[5]),y=(yyy2+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[5], sub$Q10[5], sub$Q90[5], sub$Q90[5]),y=(yyy3+3), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[6], sub$Q2.5[6], sub$Q97.5[6], sub$Q97.5[6]),y=(yyy1+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[6], sub$Q5[6], sub$Q95[6], sub$Q95[6]),y=(yyy2+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[6], sub$Q10[6], sub$Q90[6], sub$Q90[6]),y=(yyy3+4), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[7], sub$Q2.5[7], sub$Q97.5[7], sub$Q97.5[7]),y=(yyy1+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[7], sub$Q5[7], sub$Q95[7], sub$Q95[7]),y=(yyy2+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[7], sub$Q10[7], sub$Q90[7], sub$Q90[7]),y=(yyy3+5), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[8], sub$Q2.5[8], sub$Q97.5[8], sub$Q97.5[8]),y=(yyy1+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[8], sub$Q5[8], sub$Q95[8], sub$Q95[8]),y=(yyy2+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[8], sub$Q10[8], sub$Q90[8], sub$Q90[8]),y=(yyy3+6), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[9], sub$Q2.5[9], sub$Q97.5[9], sub$Q97.5[9]),y=(yyy1+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[9], sub$Q5[9], sub$Q95[9], sub$Q95[9]),y=(yyy2+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[9], sub$Q10[9], sub$Q90[9], sub$Q90[9]),y=(yyy3+7), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[10], sub$Q2.5[10], sub$Q97.5[10], sub$Q97.5[10]),y=(yyy1+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[10], sub$Q5[10], sub$Q95[10], sub$Q95[10]),y=(yyy2+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[10], sub$Q10[10], sub$Q90[10], sub$Q90[10]),y=(yyy3+8), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[11], sub$Q2.5[11], sub$Q97.5[11], sub$Q97.5[11]),y=(yyy1+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[11], sub$Q5[11], sub$Q95[11], sub$Q95[11]),y=(yyy2+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[11], sub$Q10[11], sub$Q90[11], sub$Q90[11]),y=(yyy3+9), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[14], sub$Q2.5[14], sub$Q97.5[14], sub$Q97.5[14]),y=(yyy1+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[14], sub$Q5[14], sub$Q95[14], sub$Q95[14]),y=(yyy2+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[14], sub$Q10[14], sub$Q90[14], sub$Q90[14]),y=(yyy3+10), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[13], sub$Q2.5[13], sub$Q97.5[13], sub$Q97.5[13]),y=(yyy1+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[13], sub$Q5[13], sub$Q95[13], sub$Q95[13]),y=(yyy2+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[13], sub$Q10[13], sub$Q90[13], sub$Q90[13]),y=(yyy3+11), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[15], sub$Q2.5[15], sub$Q97.5[15], sub$Q97.5[15]),y=(yyy1+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[15], sub$Q5[15], sub$Q95[15], sub$Q95[15]),y=(yyy2+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[15], sub$Q10[15], sub$Q90[15], sub$Q90[15]),y=(yyy3+12), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[12], sub$Q2.5[12], sub$Q97.5[12], sub$Q97.5[12]),y=(yyy1+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[12], sub$Q5[12], sub$Q95[12], sub$Q95[12]),y=(yyy2+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[12], sub$Q10[12], sub$Q90[12], sub$Q90[12]),y=(yyy3+13), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q2.5[16], sub$Q2.5[16], sub$Q97.5[16], sub$Q97.5[16]),y=(yyy1+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q5[16], sub$Q5[16], sub$Q95[16], sub$Q95[16]),y=(yyy2+14), col = 1, border = 1,lwd = 1)
polygon(x=c(sub$Q10[16], sub$Q10[16], sub$Q90[16], sub$Q90[16]),y=(yyy3+14), col = 1, border = 1,lwd = 1)
abline(h=4.5,lty=2)
abline(h=8.5,lty=2)
abline(h=10.5,lty=2)
abline(h=14.5,lty=2)

dev.off()
##

