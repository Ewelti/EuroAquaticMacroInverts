rm(list=ls())

library(lubridate)

#load data
d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

#choose which country for this task
TaskID <- read.csv("/data/idiv_ess/Ellen/ResponseTrends_TaskIDs.csv",as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
myCountry <- TaskID$country[which(TaskID$TaskID==task.id)]
allYrs <- subset(allYrs,country==myCountry)

#choose which response for this task
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]
allYrs$Response <- allYrs[,myResponse]

#log responses that are right skewed
if(myResponse %in% c("abundance","alien_Abund","abund_nativeSpp",
                     "EPT_Abund","insect_Abund")){
  allYrs$Response <- log10(allYrs$Response+1) 
}

#order by site site year
allYrs <- allYrs[order(allYrs$year_wMissing),]

#### two-stage models ####

### fitting directly in stan #####
library(brms)
library(rstan)

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
      
      #subset by each site
      sub <- allYrsS[allYrsS$site_id == i, ]
      
      
      #centre Year - helps model convergence to center variables for the model
      mydata$cYear <- sub$year_wMissing - median(sub$year_wMissing)
      
      #or just have as an index starting from 1
      mydata$iYear <- sub$year_wMissing - min(sub$year_wMissing)+1
      
      #centre day of year
      mydata$cday_of_year <- sub$day_of_year - median(sub$day_of_year)
      
      #if sampling occurs in more than one month include a seasonal term in the model
      maxDiffDays = max(sub$day_of_year)-min(sub$day_of_year)
      
      if(maxDiffDays < 30) {
        myformula <- bf(Response ~ cYear + cday_of_year + ar(time = iYear, p = 1, cov=FALSE))
        modelfile <- "/data/idiv_ess/Ellen/stan_code.stan"
        
      } else{
        myformula <- bf(Response ~ cYear + cday_of_year + ar(time = iYear, p = 1, cov=FALSE))
        modelfile <- "/data/idiv_ess/Ellen/stan_code_seasonal.stan"
        
      }
      
      model_data <- make_standata(myformula, data = sub, 
                                  chains = n.cores,
                                  refresh = 0)
      model_data$cYear <- sub$cYear
      model_data$cday <- sub$cday_of_year
      model_data$meanResponse <- round(median(sub$Response), 1)
      model_data$sdResponse <- max(round(mad(sub$Response), 1), 2.5)
      
      #get cores
      n.chains = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
      
      #fit model in stan
      stan_model <- stan(modelfile, 
                         data = model_data, 
                         chains = n.chains,
                         iter = 3000,
                         seed = 20)
      
      #extract model fits
      modelSummary <- summary(stan_model)$summary
      modelFits <- data.frame(estimate = modelSummary[1,"mean"],
                              sd = modelSummary[1,"sd"],
                              rhat = modelSummary[1,"Rhat"])
      trend.i <- data.frame(site = i, 
                            startyear = j,
                            modelFits)
      trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
      
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
  }} ; rm(i,j)

##save output
saveRDS(trends, file=paste0("trends__",myResponse,"__",myCountry,".RDS"))

