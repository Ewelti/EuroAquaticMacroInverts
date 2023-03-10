setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

MA <- readRDS("outputs/outputs_sensitivity/siteLevel_HTMW/HTMWsiteTrends.rds")
head(MA)
dim(MA)
unique(MA$Response)
unique(MA$Country)
unique(MA$meanYr)

MA$meanYr <- as.numeric(MA$year)+4.5
MA <- subset(MA, year <2012)

##test
df <- MA[ which(MA$Response=="abundance"), ]
dim(df)
df$fctr <- factor(df$Country)
plot(NA, xlim=range(df$meanYr), ylim=range(df$estimate))
for (x in split(df, df$siteID)) lines(x$meanYr, x$estimate, col=x$fctr[1],lwd=0.01)

###############################

#### Taxon Richness
tiff(filename = "plots/Fig3_movingWindow/HTMW_siteLevel_TaxonRich.tiff", width = 12, height = 8, units = 'in', res = 600, compression = 'lzw')

##layout
layout(mat = matrix(c(1:20), 
                        nrow = 5, 
                        ncol = 4),
       heights = c(2,2,2,2,0.5,2,2,2,2,0.5,2,2,2,2,0.5,2,2,2,2,0.5),
       widths = c(0.2,2,2,2,0.2,2,2,2,0.2,2,2,2,0.2,2,2,2,0.2,2,2,2))

par(mar=c(2,2,0.4,0.8))

#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

SR <- subset(MA, Response == "spp_richness")
head(SR)
ave_SppRich <- 27.28712314
SR$sr <-(SR$estimate/ave_SppRich)*100

#Denmark
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="a, Denmark",cex=1.5)
spl <- SR[ which(SR$Country=="Denmark"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Ireland
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="d, Ireland",cex=1.5)
spl <- SR[ which(SR$Country=="Ireland"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Norway
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="g, Norway",cex=1.5)
spl <- SR[ which(SR$Country=="Norway"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Sweden
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="j, Sweden",cex=1.5)
spl <- SR[ which(SR$Country=="Sweden"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-1,cex.lab=1.3)

#Finland
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="b, Finland",cex=1.5)
spl <- SR[ which(SR$Country=="Finland"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Latvia
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="e, Latvia",cex=1.5)
spl <- SR[ which(SR$Country=="Latvia"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Portugal
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="h, Portugal",cex=1.5)
spl <- SR[ which(SR$Country=="Portugal"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#UK
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="k, UK",cex=1.5,box.lwd = 0,box.col = "white",bg = "white")
spl <- SR[ which(SR$Country=="UK"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-1,cex.lab=1.3)

#France
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="c, France",cex=1.5)
spl <- SR[ which(SR$Country=="France"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Netherlands
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="f, Netherlands",cex=1.5)
spl <- SR[ which(SR$Country=="Netherlands"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Spain
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="i, Spain",cex=1.5)
spl <- SR[ which(SR$Country=="Spain"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-11.5,cex.lab=1.3)

dev.off()
#######################################################
###############################

#### Abundance
tiff(filename = "plots/Fig3_movingWindow/HTMW_siteLevel_Abundance.tiff", width = 12, height = 8, units = 'in', res = 600, compression = 'lzw')

##layout
layout(mat = matrix(c(1:20), 
                        nrow = 5, 
                        ncol = 4),
       heights = c(2,2,2,2,0.5,2,2,2,2,0.5,2,2,2,2,0.5,2,2,2,2,0.5),
       widths = c(0.2,2,2,2,0.2,2,2,2,0.2,2,2,2,0.2,2,2,2,0.2,2,2,2))

par(mar=c(2,2,0.4,0.8))

#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Abundance (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Abundance (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Abundance (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Abundance (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

SR <- subset(MA, Response == "abundance")
head(SR)
SR$sr <-(SR$estimate)*100

#Denmark
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="a, Denmark",cex=1.5)
spl <- SR[ which(SR$Country=="Denmark"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Ireland
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="d, Ireland",cex=1.5)
spl <- SR[ which(SR$Country=="Ireland"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Norway
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="g, Norway",cex=1.5)
spl <- SR[ which(SR$Country=="Norway"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Sweden
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="j, Sweden",cex=1.5)
spl <- SR[ which(SR$Country=="Sweden"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-1,cex.lab=1.3)

#Finland
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="b, Finland",cex=1.5)
spl <- SR[ which(SR$Country=="Finland"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Latvia
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="e, Latvia",cex=1.5)
spl <- SR[ which(SR$Country=="Latvia"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Portugal
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="h, Portugal",cex=1.5)
spl <- SR[ which(SR$Country=="Portugal"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#UK
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="k, UK",cex=1.5,box.lwd = 0,box.col = "white",bg = "white")
spl <- SR[ which(SR$Country=="UK"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-1,cex.lab=1.3)

#France
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="c, France",cex=1.5)
spl <- SR[ which(SR$Country=="France"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Netherlands
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="f, Netherlands",cex=1.5)
spl <- SR[ which(SR$Country=="Netherlands"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Spain
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="i, Spain",cex=1.5)
spl <- SR[ which(SR$Country=="Spain"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-11.5,cex.lab=1.3)

dev.off()
#######################################################

#### FRic
tiff(filename = "plots/Fig3_movingWindow/HTMW_siteLevel_FRic.tiff", width = 12, height = 8, units = 'in', res = 600, compression = 'lzw')

##layout
layout(mat = matrix(c(1:20), 
                        nrow = 5, 
                        ncol = 4),
       heights = c(2,2,2,2,0.5,2,2,2,2,0.5,2,2,2,2,0.5,2,2,2,2,0.5),
       widths = c(0.2,2,2,2,0.2,2,2,2,0.2,2,2,2,0.2,2,2,2,0.2,2,2,2))

par(mar=c(2,2,0.4,0.8))

#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Func. Richness (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Func. Richness (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Func. Richness (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Func. Richness (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

SR <- subset(MA, Response == "FRic")
head(SR)
SR$sr <-(SR$estimate)*100

#Denmark
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="a, Denmark",cex=1.5)
spl <- SR[ which(SR$Country=="Denmark"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Ireland
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="d, Ireland",cex=1.5)
spl <- SR[ which(SR$Country=="Ireland"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Norway
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="g, Norway",cex=1.5)
spl <- SR[ which(SR$Country=="Norway"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Sweden
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="j, Sweden",cex=1.5)
spl <- SR[ which(SR$Country=="Sweden"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-1,cex.lab=1.3)

#Finland
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="b, Finland",cex=1.5)
spl <- SR[ which(SR$Country=="Finland"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Latvia
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="e, Latvia",cex=1.5)
spl <- SR[ which(SR$Country=="Latvia"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Portugal
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="h, Portugal",cex=1.5)
spl <- SR[ which(SR$Country=="Portugal"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#UK
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="k, UK",cex=1.5,box.lwd = 0,box.col = "white",bg = "white")
spl <- SR[ which(SR$Country=="UK"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-1,cex.lab=1.3)

#France
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="c, France",cex=1.5)
spl <- SR[ which(SR$Country=="France"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Netherlands
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="f, Netherlands",cex=1.5)
spl <- SR[ which(SR$Country=="Netherlands"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Spain
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="i, Spain",cex=1.5)
spl <- SR[ which(SR$Country=="Spain"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-11.5,cex.lab=1.3)

dev.off()
#######################################################

#### FRed
tiff(filename = "plots/Fig3_movingWindow/HTMW_siteLevel_FRed.tiff", width = 12, height = 9, units = 'in', res = 600, compression = 'lzw')

##layout
layout(mat = matrix(c(1:20), 
                        nrow = 5, 
                        ncol = 4),
       heights = c(2,2,2,2,0.5,2,2,2,2,0.5,2,2,2,2,0.5,2,2,2,2,0.5),
       widths = c(0.2,2,2,2,0.2,2,2,2,0.2,2,2,2,0.2,2,2,2,0.2,2,2,2))

par(mar=c(2,2,0.4,0.8))

#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Func. Redundancy (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Func. Redundancy (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Func. Redundancy (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Func. Redundancy (% y"^"-1", ")")), line=-0.5,cex.lab=1.4)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

SR <- subset(MA, Response == "FRed")
ave_FRed <- 0.291145846
SR$sr <-(SR$estimate/ave_FRed)*100

#Denmark
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="a, Denmark",cex=1.5)
spl <- SR[ which(SR$Country=="Denmark"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Ireland
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="d, Ireland",cex=1.5)
spl <- SR[ which(SR$Country=="Ireland"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Norway
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="g, Norway",cex=1.5)
spl <- SR[ which(SR$Country=="Norway"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Sweden
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
title(ylab=expression(paste("Taxon richness (% y"^"-1", ")")), line=2,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="j, Sweden",cex=1.5)
spl <- SR[ which(SR$Country=="Sweden"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-1,cex.lab=1.3)

#Finland
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="b, Finland",cex=1.5)
spl <- SR[ which(SR$Country=="Finland"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Latvia
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="e, Latvia",cex=1.5)
spl <- SR[ which(SR$Country=="Latvia"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Portugal
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="h, Portugal",cex=1.5)
spl <- SR[ which(SR$Country=="Portugal"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#UK
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="k, UK",cex=1.5,box.lwd = 0,box.col = "white",bg = "white")
spl <- SR[ which(SR$Country=="UK"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-1,cex.lab=1.3)

#France
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="c, France",cex=1.5)
spl <- SR[ which(SR$Country=="France"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Netherlands
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="f, Netherlands",cex=1.5)
spl <- SR[ which(SR$Country=="Netherlands"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

#Spain
plot(NA,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=range(SR$sr), xlim=c(2004.9,2015.1))
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
legend("topright", bty="n", legend="i, Spain",cex=1.5)
spl <- SR[ which(SR$Country=="Spain"), ]
for (x in split(spl, spl$siteID)) lines(x$meanYr, x$sr, col="black",lwd=0.01)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-11.5,cex.lab=1.3)

dev.off()
#######################################################
###############################


####
#####################
###############################################
################################################