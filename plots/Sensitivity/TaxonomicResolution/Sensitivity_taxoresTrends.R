##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

# attach data
taxoresTr <- readRDS("outputs/outputs_sensitivity/TaxonomicSeason/sensitiv_taxonresTrends.RDS")
head(taxoresTr)
unique(taxoresTr$Response)

tiff(filename = "plots/Sensitivity/Sensitiv_taxoresTrends.tiff", width = 9, height = 10, units = 'in', res = 600, compression = 'lzw')

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
sub <- subset(taxoresTr, Response == "spp_richness")
yy <- c(1.5,2,2.5)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.4,0.9), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("species", "genus", "family"), las=1,cex.axis=1.3)
title(xlab="Taxon richness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### Shannon's H #####
sub <- subset(taxoresTr, Response == "shannonsH")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.04,0.025), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("species", "genus", "family"), las=1,cex.axis=1.3)
title(xlab="Shannon's H Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FRic #####
sub <- subset(taxoresTr, Response == "FRic")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.03,0.04), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("species", "genus", "family"), las=1,cex.axis=1.3)
title(xlab="Func. richness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### RaoQ #####
sub <- subset(taxoresTr, Response == "RaoQ")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.45,0.35), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("species", "genus", "family"), las=1,cex.axis=1.3)
title(xlab="Rao's Q Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### spp_richness_rarefied #####
sub <- subset(taxoresTr, Response == "spp_rich_rare")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.35,0.36), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Rarefied richness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### E10 #####
sub <- subset(taxoresTr, Response == "E10")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.015,0.01), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Shannon's evenness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FuncRedundancy #####
sub <- subset(taxoresTr, Response == "FRed")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.004,0.003), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. redundancy Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FEve #####
sub <- subset(taxoresTr, Response == "FEve")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.015,0.007), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. evenness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### abundance #####
sub <- subset(taxoresTr, Response == "abundance")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.025,0.05), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Abundance Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### turnover #####
sub <- subset(taxoresTr, Response == "turnover")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.015,0.015), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Turnover Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FDiv #####
sub <- subset(taxoresTr, Response == "FDiv")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.04,0.03), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. divergence Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### F_to #####
sub <- subset(taxoresTr, Response == "F_to")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(-0.07,0.05), cex=2,ylim=c(1.3,2.7))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. turnover Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
#legend("topleft",legend=c("family:     n = 517", "genus:    n = 537", "species:  n = 762"),bty="n",cex=2)

##
dev.off()
##
