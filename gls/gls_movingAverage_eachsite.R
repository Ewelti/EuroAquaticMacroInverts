#loop to catch each site's trends in specified intervals to calcualte moving averages

# load libraries
library(nlme)

#load data
d1 <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears.csv", header=T) # change file name according to the time series to be analyzed
allYrs <- d1[!is.na(d1$site_id_wMissing),]
allYrs$turnover<-as.numeric(allYrs$turnover) # make numeric
head(allYrs)

#centre Year - helps model convergence to center variables for the model
allYrs$cYear <- allYrs$year_wMissing - median(allYrs$year_wMissing)
summary(allYrs$cYear)
#or just have as an index starting from 1
allYrs$iYear <- allYrs$year_wMissing - min(allYrs$year_wMissing)+1
summary(allYrs$iYear)

#centre day of year
allYrs$cday_of_year <- allYrs$day_of_year - median(allYrs$day_of_year,na.rm=T)

###############################################################

##spp rich
years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(spp_richness ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(spp_richness ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "sppRich_movingAve_5yr.csv")

####

###############################################################

#spp_richness rarefied
years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(spp_rich_rare ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(spp_rich_rare ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "sppRichRare_movingAve_5yr.csv")

####
##abund
allYrs$logabun <- log10(allYrs$abundance)
years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(logabun ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(logabun ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "logabun_movingAve_5yr.csv")
##############################################################
#S_PIE
allYrs$logS_PIE <- log10(allYrs$S_PIE)

years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(logS_PIE ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(logS_PIE ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "logS_PIE_movingAve_5yr.csv")

####
##############################################################
#shannonsH

years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(shannonsH ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(shannonsH ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "shannonsH_movingAve_5yr.csv")

####
#######################################################

##turnover
years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(turnover ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(turnover ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "turnover_movingAve_5yr.csv")

####

##functional richness
allYrs$logFR <- log10(allYrs$FRic)
years <- 2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(logFR ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(logFR ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "logFuncRich_movingAve_2016.csv")

####
#######################################################

##F_to
years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(F_to ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(F_to ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "F_to_movingAve_5yr.csv")

####
#######################################################

##FEve
years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(FEve ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(FEve ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "FEve_movingAve_5yr.csv")

#################################################################
##FDiv
years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(FDiv ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(FDiv ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "FDiv_movingAve_5yr.csv")

####
#################################################################
##FDis
years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(FDis ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(FDis ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "FDis_movingAve_5yr.csv")

####
#################################################################
##RaoQ
years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      allYrsS <- subset(allYrs, year_wMissing >= j) #restrict to time period of interest
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id, #restrict to study with sufficient data in this time period - 5 years
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
      if(maxDiffDays < 30) {
        trend.i <- summary(gls(RaoQ ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      } else{
        trend.i <- summary(gls(RaoQ ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2,1]
      }
      trend.i <- data.frame(site = i, startyear = j,trend.i)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

write.csv(trends, "RaoQ_movingAve_5yr.csv")

####
