rm(list=ls())

library(lubridate)

#load data
d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

#make turnover numeric
allYrs$turnover <- as.numeric(allYrs$turnover)

#organise task IDs
TaskID <- read.csv("/data/idiv_ess/Ellen/ResponseTrends_TaskIDs.csv",as.is=T)
TaskID <- subset(TaskID, Response %in% c("alien_SppRich",
                                          "EPT_SppRich",
                                          "insect_SppRich",
                                          "insect_Abund",
                                         "EPT_Abund",
                                         "alien_Abund"))

TaskID$TaskID <- 1:nrow(TaskID)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
nrow(TaskID)#132

#choose which country for this task
myCountry <- TaskID$country[which(TaskID$TaskID==task.id)]
allYrs <- subset(allYrs,country==myCountry)

#choose which response for this task
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]
allYrs$Response <- allYrs[,myResponse]
allYrs$Response <- round(allYrs$Response)

#order by site site year
allYrs <- allYrs[order(allYrs$year_wMissing),]

#### two-stage models ####

### fitting directly in stan #####
library(brms)
library(rstan)

# write a function to consider day of year if sampling is more than 30 days apart
fitStanModel <- function(mydata){
  
  #centre Year - helps model convergence to center variables for the model
  mydata$cYear <- mydata$year_wMissing - median(mydata$year_wMissing)

  #or just have as an index starting from 1
  mydata$iYear <- mydata$year_wMissing - min(mydata$year_wMissing)+1
  
  #scale day of year
  mydata$cday_of_year <- (mydata$day_of_year - mean(mydata$day_of_year))/sd(mydata$day_of_year)
  
  #if sampling occurs in more than one month include a seasonal term in the model
  maxDiffDays = max(mydata$day_of_year)-min(mydata$day_of_year)
  
  #test zero inflation
  pois_data <-mydata$Response[!is.na(mydata$Response)]
  lambda_est <- mean(pois_data)
  p0_tilde <- exp(-lambda_est)
  p0_tilde
  n0 <- sum(1*(!(pois_data >0)))
  n <- length(pois_data)
  # number of observtions 'expected' to be zero
  n*p0_tilde
  #now lets perform the JVDB score test 
  numerator <- (n0 -n*p0_tilde)^2
  denominator <- n*p0_tilde*(1-p0_tilde) - n*lambda_est*(p0_tilde^2)
  test_stat <- numerator/denominator
  pvalue <- pchisq(test_stat,df=1, ncp=0, lower.tail=FALSE)
   
  
  if(maxDiffDays < 30 & pvalue <0.1) {
    myformula <- bf(Response ~ cYear + ar(time = iYear, p = 1, cov=FALSE))
    modelfile <- "/data/idiv_ess/Ellen/stan_code_zero.stan"
    
  } else if (maxDiffDays >= 30 & pvalue <0.1) {
    myformula <- bf(Response ~ cYear + cday_of_year + ar(time = iYear, p = 1, cov=FALSE))
    modelfile <- "/data/idiv_ess/Ellen/stan_code_seasonal_zero.stan"
    
  } else if (maxDiffDays < 30 & pvalue >0.1) {
    myformula <- bf(Response ~ cYear + cday_of_year + ar(time = iYear, p = 1, cov=FALSE))
    modelfile <- "/data/idiv_ess/Ellen/stan_code.stan"
    
  } else if (maxDiffDays >= 30 & pvalue >0.1) {
    myformula <- bf(Response ~ cYear + cday_of_year + ar(time = iYear, p = 1, cov=FALSE))
    modelfile <- "/data/idiv_ess/Ellen/stan_code_seasonal.stan"
    
  }

  # #get model code - these are saved in the above files
  # prior1 = c(set_prior("normal(0,1)", class = "b"))
  # model_code <- make_stancode(myformula, data = mydata, 
  #                              chains = 4,
  #                              family = zero_inflated_poisson(),
  #                              prior = prior1, 
  #                              refresh = 0)
  
  

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
                              family = zero_inflated_poisson(),chains = 4)
  model_data$cYear <- mydata$cYear
  model_data$cday <- mydata$cday_of_year
  model_data$meanResponse <- round(median(mydata$Response), 1)
  model_data$sdResponse <- max(round(mad(mydata$Response), 1), 2.5)
  
  #fit model in stan
  stan_model <- stan(modelfile, 
                     data = model_data, 
                     chains = 4,
                     iter = 5000,
                     init = "0",
                     control = list(adapt_delta = 0.95, 
                                    max_treedepth = 12))
  
  #extract model fits
  modelSummary <- summary(stan_model)$summary
  modelFits <- data.frame(estimate = modelSummary[1,"mean"],
                          sd = modelSummary[1,"sd"],
                          rhat = modelSummary[1,"Rhat"])
  
  }
  return(modelFits)
  
}

# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
rstan_options(auto_write = FALSE)
options(mc.cores = cpus_per_task)

#the model is called in the function below

#loop for all sites
allsites <- sort(unique(allYrs$site_id))

trends <- lapply(allsites, function(x){
  fitStanModel(subset(allYrs, site_id == x))
})

trends <- data.frame(do.call(rbind, trends))
trends$siteID <- allsites

saveRDS(trends, file=paste0("trends__",myResponse,"__",myCountry,".RDS"))
