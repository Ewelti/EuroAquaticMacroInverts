##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

###############################################
#### import site/yr summary info
sites <- read.csv("outputs/All_siteLevel_and_glmOutput.csv", header=T) # change file name according to the time series to be analyzed
CT <- data.frame(site_id=sites$site, TaxRes=sites$TaxonomicRes, country=sites$Country)
head(CT)

#### import site level estimates
ma1 <- readRDS("outputs/outputs_movingWindow/stanTrends_site_level_movingaverages.rds")
ma <- merge(ma1,CT,by="site_id", all=T) ##add in country and taxonomic resolution information
ma <- ma[complete.cases(ma), ] ## delete sites with no estimates
ma$MeanYr <- as.numeric(ma$StartYear) + 4.5
head(ma)
#### subset data into response variables
ab <- ma[ma$Response == "abundance",]
sr <- ma[ma$Response == "spp_richness",]
fred <- ma[ma$Response == "FRed",]
fric <- ma[ma$Response == "FRic",]

#######################################
#### check the proportion of sites with 
#### positive versus negative trends over time

ab_pos <- fric[fric$estimate > 0,]
ab_neg <- fric[fric$estimate < 0,]
head(ab_pos)

library(dplyr)
co_ab_pos <- ab_pos %>% count(StartYear,country)
co_ab_neg <- ab_neg %>% count(StartYear,country)
names(co_ab_pos)[names(co_ab_pos) == 'n'] <- 'n_pos'
names(co_ab_neg)[names(co_ab_neg) == 'n'] <- 'n_neg'

mo <- co_ab_pos %>% right_join(co_ab_neg, by=c("StartYear","country"))
mo[is.na(mo)] <- 0
head(mo)

mo$ab_posnegR <- mo$n_pos/(mo$n_pos+mo$n_neg)
mo$MeanYear <- as.numeric(mo$StartYear) + 4.5
max(na.omit(mo$ab_posnegR))

count_sign_mw <- mo[which(mo$MeanYear > 1994 & mo$MeanYear < 2016),]
count_sign_mw$total_sites <- (as.numeric(count_sign_mw$n_pos) + as.numeric(count_sign_mw$n_neg))
head(count_sign_mw)
co_list <-sort(unique(count_sign_mw$country))
co_list <- c("Austria",co_list)
length(co_list)

count_sign_mw$transnum <- (log10(count_sign_mw$total_sites))/2.5
count_sign_mw <- count_sign_mw[order(count_sign_mw$StartYear),] 

tiff(filename = "plots/Fig3_movingWindow/proportionPosNeg_byCountry_FuncRichness.tiff", width = 15, height = 10, units = 'in', res = 600, compression = 'lzw')

par(mfrow=c(4,5),mar=c(3.5,5,0.4,4))

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[1],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("topright", bty="n", legend="a, Austria",cex=1.5)
legend("topleft", bty="n", legend=c("trend ratio", "site count"), pch=19, col=c(1,"cyan4"))
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[2],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("bottomleft", bty="n", legend="b, Belgium",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[3],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("topright", bty="n", legend="c, Bulgaria",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[4],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("topright", bty="n", legend="d, Cyprus",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[5],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("bottomright", bty="n", legend="e, Denmark",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[6],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("topright", bty="n", legend="f, Estonia",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[7],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("bottomright", bty="n", legend="g, Finland",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[8],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("bottomright", bty="n", legend="h, France",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[9],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("bottomright", bty="n", legend="i, Germany",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[10],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("bottomright", bty="n", legend="j, Hungary",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[11],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("topright", bty="n", legend="k, Ireland",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[12],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("topright", bty="n", legend="l, Italy",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[13],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("topright", bty="n", legend="m, Latvia",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[14],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("bottomleft", bty="n", legend="n, Luxembourg",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[15],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("bottomright", bty="n", legend="o, Netherlands",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[16],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("topright", bty="n", legend="p, Norway",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[17],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("topright", bty="n", legend="q, Portugal",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[18],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("bottomright", bty="n", legend="r, Spain",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[19],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("topright", bty="n", legend="s, Sweden",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)

count_sign_mw_co <- count_sign_mw[count_sign_mw$country == co_list[20],]
plot(0,0, ylim=c(0,1),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="grey80", border = NA)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, pch=19)
points (count_sign_mw_co$ab_posnegR ~count_sign_mw_co$MeanYear, type="l")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, type="l", col="cyan4")
points (count_sign_mw_co$transnum ~ count_sign_mw_co$MeanYear, pch=19, col="cyan4")
box(lwd=2)
title(ylab="Proportion of sites with positive trends", line=2.9)
title(xlab="Mean year of moving window", line=2)
legend("bottomright", bty="n", legend="t, UK",cex=1.5)
marks <- c(2,5,10,25,50,100,200,350)
axis(4, at=c(0.120411998,0.279588002,0.4,0.559176003,0.679588002,0.8,0.920411998,1.017627218),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=3, las=3, cex=0.7)
####################
dev.off()
#########################