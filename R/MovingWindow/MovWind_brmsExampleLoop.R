####brms model for moving average

rm(list=ls())

library(brms)
library(lubridate)

#load data
d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

## only one response variable here
## e.g. spp_richness

years <- 1968:2016
trends <- NULL
for(i in unique(allYrs$site_id)){
  for(j in unique(years)) {
    tryCatch({
      
      timespan <- 10 #window length
      
      #restrict to time period of interest
      allYrsS <- subset(allYrs, year_wMissing >= j) 
      allYrsS <- subset(allYrsS, year_wMissing < (j+timespan))
      
      #restrict to study with sufficient data in this time period - 5 years
      study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id,
                              function(x)length(unique(x)))
      allYrsS <- subset(allYrsS, site_id %in% 
                          names(study_periods)[study_periods>=5])
      
      sub <- allYrsS[allYrsS$site_id == i, ]
      
      if(maxDiffDays < 30) {
        myformula <- bf(spp_richness ~ cYear + ar(time = iYear, p = 1, cov=TRUE))
      } else{
        myformula <- bf(spp_richness ~ cday_of_year + cYear + ar(time = iYear, p = 1, cov=TRUE))
      }
      
      #fit model
      fit1 <- brm(myformula, data = sub, 
                  family = gaussian(), 
                  prior = prior1, 
                  refresh = 0)
      
      #extract model fits
      modelSummary <- fixef(fit1, pars="cYear")[1, c(1:4)]
      trend.i <- data.frame(site = i, 
                            startyear = j,
                            t(modelSummary))
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
      
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

## save results
write.csv(trends, "SppRichMovingAve_brms")
