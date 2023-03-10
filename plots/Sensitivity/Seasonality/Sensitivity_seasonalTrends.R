##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

# attach data
seasTr <- readRDS("outputs/outputs_sensitivity/TaxonomicSeason/sensitiv_seasonTrends.RDS")
head(seasTr)
unique(seasTr$Response)

tiff(filename = "plots/Sensitivity/Sensitiv_seasonTrends.tiff", width = 9, height = 10, units = 'in', res = 600, compression = 'lzw')

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
sub <- subset(seasTr, Response == "spp_richness")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.5,0.55), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("spring", "summer", "fall", "winter"), las=1,cex.axis=1.3)
title(xlab="Taxon richness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### Shannon's H #####
sub <- subset(seasTr, Response == "shannonsH")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.03,0.02), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("spring", "summer", "fall", "winter"), las=1,cex.axis=1.3)
title(xlab="Shannon's H Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FRic #####
sub <- subset(seasTr, Response == "FRic")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.06,0.03), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("spring", "summer", "fall", "winter"), las=1,cex.axis=1.3)
title(xlab="Func. richness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### RaoQ #####
sub <- subset(seasTr, Response == "RaoQ")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.34,0.6), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("spring", "summer", "fall", "winter"), las=1,cex.axis=1.3)
title(xlab="Rao's Q Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### spp_richness_rarefied #####
sub <- subset(seasTr, Response == "spp_rich_rare")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.2,0.2), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Rarefied richness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### E10 #####
sub <- subset(seasTr, Response == "E10")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.01,0.01), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Shannon's evenness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FuncRedundancy #####
sub <- subset(seasTr, Response == "FRed")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.001,0.0025), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. redundancy Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FEve #####
sub <- subset(seasTr, Response == "FEve")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.005,0.015), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. evenness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### abundance #####
sub <- subset(seasTr, Response == "abundance")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.04,0.02), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Abundance Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### turnover #####
sub <- subset(seasTr, Response == "turnover")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.015,0.01), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Turnover Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FDiv #####
sub <- subset(seasTr, Response == "FDiv")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.02,0.08), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. divergence Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### F_to #####
sub <- subset(seasTr, Response == "F_to")
yy <- c(3,2,1,4)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.045,0.04), cex=2,ylim=c(0.7,4.3))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. turnover Est.", line=2.4,cex.lab=1.3)
box(lwd=2)


#plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
#legend("topleft",legend=c("winter:     n = 5", "spring:     n = 623", "summer:  n = 473", "fall:          n = 715"),bty="n",cex=2)

##
dev.off()
##
