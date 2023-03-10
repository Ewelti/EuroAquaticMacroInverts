##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

# attach data
taxoresTr <- read.csv("outputs/outputs_sensitivity/TaxonomicSeason/SplitSensitivityTrends.csv")
head(taxoresTr)
unique(taxoresTr$response)

tiff(filename = "plots/Sensitivity/Sensitiv_Split_taxoresTrends.tiff", width = 9, height = 10, units = 'in', res = 600, compression = 'lzw')

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
sub <- subset(taxoresTr, response == "SppRich")
yy <- c(1.5,2,2.5,3)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("species", "genus", "family", "all"), las=1,cex.axis=1.3)
title(xlab="Taxon richness Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#### Shannon's H #####
sub <- subset(taxoresTr, response == "shannonsH")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("species", "genus", "family", "all"), las=1,cex.axis=1.3)
title(xlab="Shannon's H Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FRic #####
sub <- subset(taxoresTr, response == "FRic")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("species", "genus", "family", "all"), las=1,cex.axis=1.3)
title(xlab="Functional richness Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#### RaoQ #####
sub <- subset(taxoresTr, response == "RaoQ")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("species", "genus", "family", "all"), las=1,cex.axis=1.3)
title(xlab="Rao's Q Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#### spp_richness_rarefied #####
sub <- subset(taxoresTr, response == "SppRichRare")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Rarefied richness Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#### E10 #####
sub <- subset(taxoresTr, response == "E10")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Shannon's evenness Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FuncRedundancy #####
sub <- subset(taxoresTr, response == "FRed")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Functional redundancy Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FEve #####
sub <- subset(taxoresTr, response == "FEve")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Functional evenness Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#### abundance #####
sub <- subset(taxoresTr, response == "abundance")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Abundance Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#### turnover #####
sub <- subset(taxoresTr, response == "turnover")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Turnover Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FDiv #####
sub <- subset(taxoresTr, response == "FDiv")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Functional divergence Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#### F_to #####
sub <- subset(taxoresTr, response == "Fto")
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(1.3,3.2))
polygon(x=c(-100,-100,0,0),
        y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Functional turnover Estimate", line=2.4,cex.lab=1.3)
box(lwd=2)

#plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
#legend("topleft",legend=c("family:     n = 517", "genus:    n = 537", "species:  n = 762"),bty="n",cex=2)

##
dev.off()
##
