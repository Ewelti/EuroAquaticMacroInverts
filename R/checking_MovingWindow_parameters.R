##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

###############################################
#### import site/yr summary info
sites <- read.csv("/outputs/All_siteLevel_and_glmOutput.csv", header=T) # change file name according to the time series to be analyzed
CT <- data.frame(site_id=sites$site, TaxRes=sites$TaxonomicRes, country=sites$Country)
head(CT)

#### import year estimates
ma_ests <- read.csv("/outputs/movingAve_YrEsts.csv")
#### subset data into response variables
ab_est <- ma_ests[ma_ests$Response == "abundance",]
sr_est <- ma_ests[ma_ests$Response == "spp_richness",]

#### import site level estimates
ma1 <- readRDS("/outputs/stanTrends_site_level_movingaverages.rds")
ma1 <- merge(ma1,CT,by="site_id", all=T) ##add in country and taxonomic resolution information
ma <- ma[complete.cases(ma), ] ## delete sites with no estimates
ma$MeanYr <- as.numeric(ma$StartYear) + 4.5
head(ma_ests)
head(ma)
#### subset data into response variables
ab <- ma[ma$Response == "abundance",]
sr <- ma[ma$Response == "spp_richness",]

#install.packages("data.table")
library(data.table)

#################################################
##### check changes in number of sites per country and per taxonomic resolution over time

# number of sites per county per year
ab_or <- ab[order(ab$MeanYr),]
DT <- data.table(ab_or)
CountryCount <- DT[, .(site_num = length(unique(site_id))), by = list(MeanYr,country)]
head(CountryCount)
unique(CountryCount$country)

# number of sites per taxonomic resolution per year
TaxResCount <- DT[, .(site_num = length(unique(site_id))), by = list(MeanYr,TaxRes)]
head(TaxResCount)

#line plots
cc <- CountryCount[ which(CountryCount$MeanYr > 1994 & CountryCount$MeanYr < 2016), ]
cc$fctr <- factor(cc$county)
plot(NA, ylim=range(cc$site_num),xlim= c(1994,2016))
for (x in split(cc, cc$country)) lines(x$MeanYr, x$site_num, col=x$fctr[1],lwd=0.01)

TaxResCount <- TaxResCount[order(TaxResCount$MeanYr),]
tr <- TaxResCount[ which(TaxResCount$MeanYr > 1994 & TaxResCount$MeanYr < 2016), ]
tr$fctr <- factor(tr$TaxRes)
plot(NA, ylim=range(tr$site_num),xlim= c(1994,2016))
for (x in split(tr, tr$TaxRes)) lines(x$MeanYr, x$site_num, col=x$fctr[1],lwd=0.01)

##barplots
library(tidyr)
cut <- CountryCount[ which(CountryCount$MeanYr > 1994 & CountryCount$MeanYr < 2016), ]
ccw <- spread(cut, MeanYr, site_num)
ccw [is.na(ccw )] <- 0
ccw

oth <- ccw[1,2:23] + ccw[3,2:23] + ccw[4,2:23] + ccw[6,2:23] + ccw[7,2:23] + ccw[11,2:23] + ccw[12,2:23] + ccw[13,2:23] + ccw[14,2:23] + ccw[17,2:23]
lab <- data.frame(country="Other",oth)
names(lab) <- sub("^X", "", names(lab))

ccww <- rbind(lab,ccw[20,],ccw[19,],ccw[18,],ccw[16,],ccw[15,],ccw[10,],ccw[9,],ccw[8,],ccw[5,],ccw[2,])
ccww[2,1] = "UK"

# # bargraph # #
tiff(filename = "plots/MA_countriesPerYr.tiff", width =5, height = 8, units = 'in', res = 600, compression = 'lzw')

par(mgp=c(3,0,-1),mar=c(4,4,0.4,0.2)+0.1)

colors <- c("black","lightslateblue","lightskyblue", "darkblue","darkseagreen2", "forestgreen",
"brown1","darkorange","red3","peachpuff", "saddlebrown")

barplot(as.matrix(ccww), las=2, col= colors, border=colors, xaxt='n')

place <- seq(1.6,26.8, length = 22)
marks <- c(1994:2015)
axis(side = 1, at=place, labels=marks,las=2,tick=F,line=F)
legend(x=2,y=1400, rev(rownames(sites)), fill = rev(colors), border = rev(colors), bty = "n",cex=1.4)
title(ylab="Number of sites", line=2.2,cex.lab=1.5)
title(xlab="Mean year of moving window", line=2.2,cex.lab=1.5)

##
dev.off()
##


##########################################
## check changes in estimate error over time

ab_est$MeanYear <- as.numeric(ab_est$StartYear) + 4.5
ab_est_mw <- ab_est[ab_est$MeanYear > 1994,]
ab_est_mw <- ab_est_mw[ab_est_mw$MeanYear < 2016,]
head(ab_est_mw)

sr_est$MeanYear <- as.numeric(sr_est$StartYear) + 4.5
sr_est_mw <- sr_est[sr_est$MeanYear > 1994,]
sr_est_mw <- sr_est_mw[sr_est_mw$MeanYear < 2016,]

par(mfrow=c(1,2))
max(sr_est_mw$Est.Error)
plot(sr_est_mw$Est.Error ~sr_est_mw$MeanYear, ylim=c(0,0.35),xlim= c(1994,2016), pch=19,
ylab="Estimated error of taxon richness trend", xlab="Mean year of moving window")
points (sr_est_mw$Est.Error ~sr_est_mw$MeanYear, type="l")

max(ab_est_mw$Est.Error)
plot(ab_est_mw$Est.Error ~ab_est_mw$MeanYear, ylim=c(0,0.025),xlim= c(1994,2016), pch=19,
ylab="Estimated error of abundance trend", xlab="Mean year of moving window")
points (ab_est_mw$Est.Error ~ab_est_mw$MeanYear, type="l")

#######################################
#### check the proportion of sites with 
#### positive versus negative trends over time

ab_pos <- ab[ab$estimate > 0,]
ab_neg <- ab[ab$estimate < 0,]

sr_pos <- sr[sr$estimate > 0,]
sr_neg <- sr[ab$estimate < 0,]

DT <- data.table(ab_pos)
sitecount_ab_pos <- DT[, .(site_num_Abund_pos = length(unique(site_id))), by = StartYear]
DT <- data.table(ab_neg)
sitecount_ab_neg <- DT[, .(site_num_Abund_neg = length(unique(site_id))), by = StartYear]
sitecount_ab_neg

DT_sr <- data.table(sr_pos)
sitecount_sr_pos <- DT_sr[, .(site_num_SRich_pos = length(unique(site_id))), by = StartYear]
sitecount_sr_pos
DT_sr <- data.table(sr_neg)
sitecount_sr_neg <- DT_sr[, .(site_num_SRich_neg = length(unique(site_id))), by = StartYear]
sitecount_sr_neg

MoAv1 <- merge(sitecount_ab_pos,sitecount_ab_neg,by="StartYear", all=T)
MoAv2 <- merge(MoAv1,sitecount_sr_pos,by="StartYear", all=T)
count_sign <- merge(MoAv2,sitecount_sr_neg,by="StartYear", all=T)
head(count_sign)

count_sign$ab_posnegR <- count_sign$site_num_Abund_pos/count_sign$site_num_Abund_neg
count_sign$sr_posnegR <- count_sign$site_num_SRich_pos/count_sign$site_num_SRich_neg
count_sign$MeanYear <- as.numeric(count_sign$StartYear) + 4.5

count_sign_mw <- count_sign[count_sign$MeanYear > 1994,]
count_sign_mw <- count_sign_mw[count_sign_mw$MeanYear < 2016,]
count_sign_mw$total_sites <- (as.numeric(count_sign_mw$site_num_Abund_pos) + as.numeric(count_sign_mw$site_num_Abund_neg))

transnum <- log10(count_sign_mw$total_sites)-1.6

par(mfrow=c(1,2),mar=c(3,5,0.2,4))
plot(0,0, ylim=c(0.7,1.7),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 1, 1, -100), col ="grey80", border = NA)
points (transnum ~ count_sign_mw$MeanYear, type="l", col="cyan4")
points (transnum ~ count_sign_mw$MeanYear, pch=19, col="cyan4")
points (count_sign_mw$sr_posnegR ~count_sign_mw$MeanYear, pch=19)
points (count_sign_mw$sr_posnegR ~count_sign_mw$MeanYear, type="l")
box(lwd=2)
title(ylab="# sites with positive/negative", line=3.2)
title(ylab="taxon richness trends", line=2.2)
title(xlab="Mean year of moving window", line=2)
marks <- c(250, 500, 750, 1000, 1250, 1500)
axis(4, at=c(0.797940009, 1.098970004, 1.275061263, 1.4, 1.496910013, 1.576091259),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=2.5, las=3)

plot(0,0, ylim=c(0.7,1.7),xlim= c(1994,2016), ylab="", xlab="", las=1)
polygon(x = c(1, 1, 2040, 2040), y = c(-100, 1, 1, -100), col ="grey80", border = NA)
points (transnum ~ count_sign_mw$MeanYear, type="l", col="cyan4")
points (transnum ~ count_sign_mw$MeanYear, pch=19, col="cyan4")
points (count_sign_mw$ab_posnegR ~count_sign_mw$MeanYear, pch=19)
points (count_sign_mw$ab_posnegR ~count_sign_mw$MeanYear, type="l")
box(lwd=2)
title(ylab="# sites with positive/negative", line=3.2)
title(ylab="abundance trends", line=2.2)
title(xlab="Mean year of moving window", line=2)
legend("topleft", bty="n", legend=c("trend ratio", "site count"), pch=19, col=c(1,"cyan4"))
axis(4, at=c(0.797940009, 1.098970004, 1.275061263, 1.4, 1.496910013, 1.576091259),cex.axis=1,labels=marks,las=1,col="cyan4")
mtext(side=4, "Number of sites", line=2.5)
###############################################################################

#### check average number of years sampled per time series
#### in each year of the moving window analysis

## attach data
DATA1_list <- read.csv("All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T)
allYrs <- DATA1_list[!is.na(DATA1_list$site_id_wMissing),]
attach(allYrs)
head(allYrs)

#### loop to count years in moving window analysis

#######WARNING: takes a few hours to run

years <- 1968:2016
count <- NULL
for(i in unique(allYrs$site_id)){
for(j in unique(years)) {
  tryCatch({
	
	timespan <- 10 #window length

      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))

 	study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to studies with sufficient data in this time period - 6 years
                          function(x)length(unique(x)))
  	allYrsS <- subset(allYrsS, site_id %in% 
                      names(study_periods)[study_periods>=6])

	sub <- allYrsS[allYrsS$site_id == i, ]
	count.i <- nrow(sub) 

    	count.i <- data.frame(site = i, startyear = j,count.i)
    	count <- rbind(count, count.i) ; rm(count.i, sub)

  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
}} ; rm(i,j)

####
count_nz <- count[count$count.i > 0,]
head(count_nz)

write.csv(count_nz, "yearcount_sitelevel_MovingWindow.csv")
count_nz <- read.csv("yearcount_sitelevel_MovingWindow.csv")


###########################################################
count_nz$myr <- as.numeric(count_nz$startyear) + 4.5
count_nz_cut <- count_nz[count_nz$startyear > 1994,]

boxplot(count_nz_cut$count.i ~count_nz_cut$myr, ylab="Sampling Years", xlab="Mean year of moving window")



####
## get mean number of years sampled per site/time series 
## with each moving window period

mean_yr_count <- aggregate(count_nz$count.i, list(count_nz$startyear), FUN = mean)
colnames(mean_yr_count) <- c("startyear", "count.i")
head(mean_yr_count)

## plot # mean years sampled over moving window years

mean_yr_count$MeanYear <- as.numeric(mean_yr_count$startyear) + 4.5
mean_yr_count_mw <- mean_yr_count[ which(mean_yr_count$MeanYear >1994 & mean_yr_count$MeanYear <2016), ]

plot(mean_yr_count_mw$count.i ~mean_yr_count_mw$MeanYear, ylim=c(6,9),xlim= c(1994,2016), pch=19,
ylab="Mean number of sampled years", xlab="Mean year of moving window")
points (mean_yr_count_mw$count.i ~mean_yr_count_mw$MeanYear, type="l")

#########################################



