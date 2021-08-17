rm(list=ls())

library(lubridate)

#load data
d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

#choose which country for this task
TaskID <- read.csv("/data/idiv_ess/Ellen/ResponseTrends_TaskIDs_rerun.csv",as.is=T)
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

#centre Year - helps model convergence to center variables for the model
allYrs$cYear <- allYrs$year_wMissing - median(allYrs$year_wMissing)
summary(allYrs$cYear)
#or just have as an index starting from 1
allYrs$iYear <- allYrs$year_wMissing - min(allYrs$year_wMissing)+1
summary(allYrs$iYear)

#centre day of year
allYrs$cday_of_year <- allYrs$day_of_year - median(allYrs$day_of_year,na.rm=T)

#### two-stage models ####

### fitting directly in stan #####
library(brms)
library(rstan)

# write a function to consider day of year if sampling is more than 30 days apart
fitStanModel <- function(mydata){
  
  #if sampling occurs in more than one month include a seasonal term in the model
  maxDiffDays = max(mydata$cday_of_year)-min(mydata$cday_of_year)
  
  if(maxDiffDays < 30) {
    myformula <- bf(Response ~ cYear + cday_of_year + ar(time = iYear, p = 1, cov=FALSE))
    
    modelfile <- "/data/idiv_ess/Ellen/stan_code.stan"
    
  } else{
    myformula <- bf(Response ~ cYear + cday_of_year + ar(time = iYear, p = 1, cov=FALSE))
    modelfile <- "/data/idiv_ess/Ellen/stan_code_seasonal.stan"
    
  }
  
#responses: spp_richness

  #get model code - these are saved in the above files
  # prior1 = c(set_prior("normal(0,10)", class = "b"))
  # model_code <- make_stancode(myformula, data = mydata, 
  #                             chains = n.cores,
  #                             prior = prior1, 
  #                             refresh = 0)
  
  

  #there is missing data, dont run the model
  if(all(is.na(mydata$Response))){
    
    modelFits <- data.frame(estimate = NA,
                            sd = NA,
                            rhat = NA)
  }else{
    
  #get model data
  
  #remove any missing values
  mydata <- subset(mydata, !is.na(Response))
      
  model_data <- make_standata(myformula, data = mydata, 
                                chains = n.cores,
                                refresh = 0)
  model_data$cYear <- mydata$cYear
  model_data$cday <- mydata$cday_of_year
    
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
  
  }
  return(modelFits)
  
}

#get cores
n.chains = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

#loop for all sites
allsites <- sort(unique(allYrs$site_id))

trends <- lapply(allsites, function(x){
  fitStanModel(subset(allYrs, site_id == x))
})

trends <- data.frame(do.call(rbind, trends))
trends$siteID <- allsites

saveRDS(trends, file=paste0("trends__",myResponse,"__",myCountry,".RDS"))

