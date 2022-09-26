##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

# attach data
JKTr <- read.csv("outputs/JackknifeTrends.csv")
head(JKTr)
unique(JKTr$Response)

tiff(filename = "plots/Sensitivity/JacknifeTrends.tiff", width = 9, height = 10, units = 'in', res = 600, compression = 'lzw')

layout(mat = matrix(c(1:16), 
                        nrow = 4, 
                        ncol = 4),
       heights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
       widths = c(0.6,2,2,2,0.6,2,2,2,0.6,2,2,2,0.6,2,2,2))
layout.show(16)

par(mar=c(4,0.4,0.4,0.4))
#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
#### Spp Richness #####
sub <- subset(JKTr, Response == "spp_richness")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=(sub$country), las=1,cex.axis=0.7)
title(xlab=expression(paste("Taxon richness ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#### Shannon's H #####
sub <- subset(JKTr, Response == "shannonsH")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=(sub$country), las=1,cex.axis=0.7)
title(xlab=expression(paste("Shannon's H ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#### FRic #####
sub <- subset(JKTr, Response == "FRic")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=(sub$country), las=1,cex.axis=0.7)
title(xlab=expression(paste("Func. richness ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#### RaoQ #####
sub <- subset(JKTr, Response == "RaoQ")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=(sub$country), las=1,cex.axis=0.7)
title(xlab=expression(paste("Rao's Q ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#### spp_richness_rarefied #####
sub <- subset(JKTr, Response == "spp_rich_rare")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
#axis(2, at=yy, labels=(sub$country), las=1,cex.axis=1.3)
title(xlab=expression(paste("Rarefied richness ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#### E10 #####
sub <- subset(JKTr, Response == "E10")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
#axis(2, at=yy, labels=(sub$country), las=1,cex.axis=1.3)
title(xlab=expression(paste("Shannon's evenness ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#### FuncRedundancy #####
sub <- subset(JKTr, Response == "FRed")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
#axis(2, at=yy, labels=(sub$country), las=1,cex.axis=1.3)
title(xlab=expression(paste("Func. redundancy ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#### FEve #####
sub <- subset(JKTr, Response == "FEve")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
#axis(2, at=yy, labels=(sub$country), las=1,cex.axis=1.3)
title(xlab=expression(paste("Func. evenness ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#### abundance #####
sub <- subset(JKTr, Response == "abundance")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
#axis(2, at=yy, labels=(sub$country), las=1,cex.axis=1.3)
title(xlab=expression(paste("Abundance ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#### turnover #####
sub <- subset(JKTr, Response == "turnover")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
#axis(2, at=yy, labels=(sub$country), las=1,cex.axis=1.3)
title(xlab=expression(paste("Turnover ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#### FDiv #####
sub <- subset(JKTr, Response == "FDiv")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
#axis(2, at=yy, labels=(sub$country), las=1,cex.axis=1.3)
title(xlab=expression(paste("Func. divergence ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#### F_to #####
sub <- subset(JKTr, Response == "F_to")
yy <- c(22:1)
plot(yy ~ sub$percEstimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$percQ2.5),max(sub$percQ97.5)), cex=2,ylim=c(0.5,22.5))
polygon(x=c(-100,-100,0,0),
        y=c(-4,29,29,-4), col = "grey80", border = "grey80")
points(yy ~ sub$percEstimate, pch=19, cex=1)
arrows(sub$percQ2.5, yy, sub$percQ97.5, yy, length=0.05, angle=90, code=3, lwd=2)
#axis(2, at=yy, labels=(sub$country), las=1,cex.axis=1.3)
title(xlab=expression(paste("Func. turnover ",y^-1,, sep = "")), line=2.4,cex.lab=1.3)
box(lwd=2)

#plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
#legend("topleft",legend=c("family:     n = 517", "genus:    n = 537", "species:  n = 762"),bty="n",cex=2)

##
dev.off()
##
