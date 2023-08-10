##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

#add pcks
library(tidyverse)
library(janitor)

# attach data
MA <- read.csv("outputs/outputs_movingWindow/movingAve_YrEsts.csv")
MA$meanYr <- (MA$StartYear+ 4.5)
head(MA)
max(MA$StartYear)
MA = subset(MA, select = -c(X) )

metricdata <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv")
head(metricdata)

#############################################
movav1 <- readRDS("outputs/outputs_movingWindow/proportions_modelSummaries.rds")
movav1 <- data.frame(movav1)
movav1$MeanYear <- (movav1$start_year+ 4.5)
colnames(movav1)[2] ="StartYear"
head(movav1)
max(movav1$StartYear)
movav1 <- movav1[ which(movav1$MeanYear < 2016), ]

##################################################

#tiff(filename = "plots/Fig3_movingWindow/Fig3_MovingAveragesProportions.tiff", width = 11, height = 6.5, units = 'in', res = 600, compression = 'lzw')
pdf("plots/Fig3_movingWindow/Fig3_MovingAveragesProportions.pdf", width = 11, height = 6.5)


##layout
layout(mat = matrix(c(1:15), 
                        nrow = 3, 
                        ncol = 5),
       heights = c(2,2,0.4),
       widths = c(0.5,2,2,2,2))

par(mar=c(1.8,2,0.4,1))
#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=expression(paste("Mean site trend (% change y"^"-1",")")), line=-2,cex.lab=1.6)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab="Proportion sites with", line=-0.6,cex.lab=1.6)
title(ylab="positive trends", line=-2,cex.lab=1.6)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

#plot for spp richness
SR <- subset(MA, Response == "spp_richness")
ave_SppRich <- mean(metricdata$spp_richness, na.rm=T)
st <-(SR[1:(nrow(SR)),3:12]/ave_SppRich)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
#SRs <- sr[ which(sr$site_num >=100), ]
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ] # this is min 7 countries
#SRs <- sr[ which(sr$StartYear >=1980), ] #this is min 5 countries

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
title(ylab=expression(paste("% y"^"-1")), line=2.4,cex.lab=1.4)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)#pch="-",
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("bottomright", bty="n", legend="Taxon richness",cex=1.5)
legend("topleft", bty="n", legend=as.expression(bquote(bold("a"))),cex=1.5)

count_sign_mw <- movav1[which(movav1$response=='spp_richness'),]
plot(0,0, ylim=c(0.3,0.8),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower90Prop,rev(count_sign_mw$upper90Prop)),col="#00000020", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower80Prop,rev(count_sign_mw$upper80Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
legend("bottomright", bty="n", legend="Taxon richness",cex=1.5)
legend("topleft", bty="n", legend=as.expression(bquote(bold("e"))),cex=1.5)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-2,cex.lab=1.3)

#plot for abundance
SR <- subset(MA, Response == "abundance")
st <-(10^(SR[1:(nrow(SR)),3:12])-1)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)#pch="-",
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("bottomright", bty="n", legend="Abundance",cex=1.5)
legend("topleft", bty="n", legend=as.expression(bquote(bold("b"))),cex=1.5)

count_sign_mw <- movav1[which(movav1$response=='abundance'),]
plot(0,0, ylim=c(0.3,0.8),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower90Prop,rev(count_sign_mw$upper90Prop)),col="#00000020", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower80Prop,rev(count_sign_mw$upper80Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
title(xlab="Mean year of moving window", line=2.5,cex.lab=1.1)
legend("bottomright", bty="n", legend="Abundance",cex=1.5)
legend("topleft", bty="n", legend=as.expression(bquote(bold("f"))),cex=1.5)


plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-2,cex.lab=1.3)

#plot for FRic
SR <- subset(MA, Response == "FRic")
st <-(10^(SR[1:(nrow(SR)),3:12])-1)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("bottomright", bty="n", legend="Functional richness",cex=1.5)
legend("topleft", bty="n", legend=as.expression(bquote(bold("c"))),cex=1.5)

count_sign_mw <- movav1[which(movav1$response=='FRic'),]
plot(0,0, ylim=c(0.3,0.8),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower90Prop,rev(count_sign_mw$upper90Prop)),col="#00000020", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower80Prop,rev(count_sign_mw$upper80Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
title(xlab="Mean year of moving window", line=2.5,cex.lab=1.1)
legend("bottomright", bty="n", legend="Functional richness",cex=1.5)
legend("topleft", bty="n", legend=as.expression(bquote(bold("g"))),cex=1.5)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-2,cex.lab=1.3)

#plot for FRed
SR <- subset(MA, Response == "FRed")
ave_FRed <- mean(metricdata$FRed, na.rm=T)
st <-(SR[1:(nrow(SR)),3:12]/ave_FRed)*100
sr <- cbind(SR$StartYear, SR$site_num, SR$meanYr, st)
names(sr)[names(sr) == 'SR$StartYear'] <- 'StartYear'
names(sr)[names(sr) == 'SR$site_num'] <- 'site_num'
names(sr)[names(sr) == 'SR$meanYr'] <- 'meanYr'
#select yrs with enough sites to be representative
SRs <- sr[ which(sr$StartYear >1989 & sr$StartYear <2012), ]

plot(SRs$Estimate~SRs$meanYr,ylab="", cex=1.5, xlab="", type="n", las=1, ylim=c(min(SRs$Q2.5),max(SRs$Q97.5)), xlim=c(min(SRs$meanYr),max(SRs$meanYr)))
title(xlab="Mean year of moving window", line=2.4,cex.lab=1.3)
polygon(x = c(0, 0, 2040, 2040), y = c(-100, 0, 0, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0, 0, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q2.5,rev(SRs$Q97.5)),col="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q5,rev(SRs$Q95)), col ="#00000020", border = NA)
polygon(x = c(SRs$meanYr,rev(SRs$meanYr)), y = c(SRs$Q10,rev(SRs$Q90)), col ="#00000020", border = NA)
points(SRs$Estimate~SRs$meanYr, pch=21,bg=1,cex=1.8)
points(SRs$Estimate~SRs$meanYr,type="l",lwd=2)
legend("bottomright", bty="n", legend="Functional redundancy",cex=1.5)
legend("topleft", bty="n", legend=as.expression(bquote(bold("d"))),cex=1.5)

count_sign_mw <- movav1[which(movav1$response=='FRed'),]
plot(0,0, ylim=c(0.3,0.8),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower90Prop,rev(count_sign_mw$upper90Prop)),col="#00000020", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower80Prop,rev(count_sign_mw$upper80Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
title(xlab="Mean year of moving window", line=2.5,cex.lab=1.1)
legend("bottomright", bty="n", legend="Functional redundancy",cex=1.5)
legend("topleft", bty="n", legend=as.expression(bquote(bold("h"))),cex=1.5)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Mean year of moving window", line=-2,cex.lab=1.3)

##

dev.off()
###############################################################################


