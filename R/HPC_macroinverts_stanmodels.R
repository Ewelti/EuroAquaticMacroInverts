library(brms)
library(lubridate)

#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

#load data
d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

#centre Year - helps model convergence to center variables for the model
allYrs$cYear <- allYrs$year_wMissing - median(allYrs$year_wMissing)
summary(allYrs$cYear)
#or just have as an index starting from 1
allYrs$iYear <- allYrs$year_wMissing - min(allYrs$year_wMissing)+1
summary(allYrs$iYear)

#centre day of year
allYrs$cday_of_year <- allYrs$day_of_year - median(allYrs$day_of_year,na.rm=T)

#### two-stage models ####

#set priors now
prior1 = c(set_prior("normal(0,0.5)", class = "ar"),
           set_prior("normal(0,10)", class = "b"))


# write a function to consider day of year if sampling is more than 30 days apart
fitStanModel <- function(mydata){
  
  #if sampling occurs in more than one month include a seasonal term in the model
  maxDiffDays = max(mydata$cday_of_year)-min(mydata$cday_of_year)
  
  if(maxDiffDays < 30) {
    myformula <- bf(spp_richness ~ cYear + ar(time = iYear, p = 1))
  } else{
    myformula <- bf(spp_richness ~ cday_of_year + cYear + ar(time = iYear, p = 1))
  }
  
  #fit model
  fit1 <- brm(myformula, data = mydata, family = poisson(), prior = prior1, refresh = 0)
  
  #extract model fits
  modelSummary <- fixef(fit1, pars="cYear")[1, c(1,2)]
  return(modelSummary)
  
}

#loop for all sites ####has problems with crashing-- trying tryCatch
trends <- NULL
for(i in unique(allYrs$site_id)){
  tryCatch({
    sub <- allYrs[allYrs$site_id == i, ]
    trend.i <- fitStanModel(sub)
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})    
} ; rm(i)
saveRDS(trends, file="trends.RDS")
