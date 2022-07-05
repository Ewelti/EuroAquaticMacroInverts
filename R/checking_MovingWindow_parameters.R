#### set working directory
setwd("C:/Users/Ellen/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs")

###############################################
#### import data and library
ma_ests <- read.csv("movingAve_YrEsts.csv")
head(ma_ests)
#### subset data into response variables
ab_est <- ma_ests[ma_ests$Response == "abundance",]
sr_est <- ma_ests[ma_ests$Response == "spp_richness",]

ma <- readRDS("stanTrends_site_level_movingaverages.rds")
head(ma)
#### subset data into response variables
ab <- ma[ma$Response == "abundance",]
sr <- ma[ma$Response == "spp_richness",]

library(data.table)

##########################################
## check changes in estimate error over time

ab_est$MeanYear <- as.numeric(ab_est$StartYear) + 4.5
ab_est_mw <- ab_est[ab_est$MeanYear > 1994,]
head(ab_est_mw)

sr_est$MeanYear <- as.numeric(sr_est$StartYear) + 4.5
sr_est_mw <- sr_est[sr_est$MeanYear > 1994,]

par(mfrow=c(1,2))

max(ab_est_mw$Est.Error)
plot(ab_est_mw$Est.Error ~ab_est_mw$MeanYear, ylim=c(0,0.025),xlim= c(1994,2019), pch=19,
ylab="Estimated error of abundance trend", xlab="Mean year of moving window")
points (ab_est_mw$Est.Error ~ab_est_mw$MeanYear, type="l")

max(sr_est_mw$Est.Error)
plot(sr_est_mw$Est.Error ~sr_est_mw$MeanYear, ylim=c(0,0.35),xlim= c(1994,2019), pch=19,
ylab="Estimated error of taxon richness trend", xlab="Mean year of moving window")
points (sr_est_mw$Est.Error ~sr_est_mw$MeanYear, type="l")

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

par(mfrow=c(1,2))
plot(count_sign_mw$sr_posnegR ~count_sign_mw$MeanYear, ylim=c(0.7,1.7),xlim= c(1994,2019), pch=19,
ylab="# sites with positive/negative trends in Taxon richness ", xlab="Mean year of moving window")
points (count_sign_mw$sr_posnegR ~count_sign_mw$MeanYear, type="l")
abline(h=1, lty=2, col="grey60")

plot(count_sign_mw$ab_posnegR ~count_sign_mw$MeanYear, ylim=c(0.7,1.7),xlim= c(1994,2019), pch=19,
ylab="# sites with positive/negative trends in Abundance", xlab="Mean year of moving window")
points (count_sign_mw$ab_posnegR ~count_sign_mw$MeanYear, type="l")
abline(h=1, lty=2, col="grey60")
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

####
## get mean number of years sampled per site/time series 
## with each moving window period

mean_yr_count <- aggregate(count_nz$count.i, list(count_nz$startyear), FUN = mean)
colnames(mean_yr_count) <- c("startyear", "count.i")
head(mean_yr_count)

## plot # mean years sampled over moving window years

mean_yr_count$MeanYear <- as.numeric(mean_yr_count$startyear) + 4.5
mean_yr_count_mw <- mean_yr_count[mean_yr_count$MeanYear > 1994,]

plot(mean_yr_count_mw$count.i ~mean_yr_count_mw$MeanYear, ylim=c(6,9),xlim= c(1994,2019), pch=19,
ylab="Mean number of sampled years", xlab="Mean year of moving window")
points (mean_yr_count_mw$count.i ~mean_yr_count_mw$MeanYear, type="l")

#########################################



