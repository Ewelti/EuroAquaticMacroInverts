##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

# attach data
decTr <- read.csv("outputs/DecadeTrends.csv")
head(decTr)
unique(decTr$Response)

tiff(filename = "plots/DecadeTrends.tiff", width = 9, height = 10, units = 'in', res = 600, compression = 'lzw')

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
sub <- subset(decTr, Response == "spp_richness")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("2000s","2010s"), las=1,cex.axis=1.3)
title(xlab="Taxon richness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### Shannon's H #####
sub <- subset(decTr, Response == "shannonsH")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("2000s","2010s"), las=1,cex.axis=1.3)
title(xlab="Shannon's H Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FRic #####
sub <- subset(decTr, Response == "FRic")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("2000s","2010s"), las=1,cex.axis=1.3)
title(xlab="Func. richness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### RaoQ #####
sub <- subset(decTr, Response == "RaoQ")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
axis(2, at=yy, labels=c("2000s","2010s"), las=1,cex.axis=1.3)
title(xlab="Rao's Q Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### spp_richness_rarefied #####
sub <- subset(decTr, Response == "spp_rich_rare")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Rarefied richness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### E10 #####
sub <- subset(decTr, Response == "E10")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Shannon's evenness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FuncRedundancy #####
sub <- subset(decTr, Response == "FRed")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. redundancy Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FEve #####
sub <- subset(decTr, Response == "FEve")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. evenness Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### abundance #####
sub <- subset(decTr, Response == "abundance")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Abundance Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### turnover #####
sub <- subset(decTr, Response == "turnover")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Turnover Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### FDiv #####
sub <- subset(decTr, Response == "FDiv")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. divergence Est.", line=2.4,cex.lab=1.3)
box(lwd=2)

#### F_to #####
sub <- subset(decTr, Response == "F_to")
sub <- sub[3:4,]
yy <- c(2,1)
plot(yy ~ sub$Estimate, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(min(sub$Q2.5),max(sub$Q97.5)), cex=2,ylim=c(0.7,2.3))
polygon(x=c(-100,-100,0,0), y=c(-4,22,22,-4), col = "grey80", border = "grey80")
points(yy ~ sub$Estimate, pch=19,cex=2)
arrows(sub$Q2.5, yy, sub$Q97.5, yy, length=0.05, angle=90, code=3, lwd=2)
title(xlab="Func. turnover Est.", line=2.4,cex.lab=1.3)
box(lwd=2)


#plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
#legend("topleft",legend=c("winter:     n = 5", "spring:     n = 623", "summer:  n = 473", "fall:          n = 715"),bty="n",cex=2)

##
dev.off()
##
