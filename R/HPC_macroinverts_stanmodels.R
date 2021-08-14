rm(list=ls())

library(lubridate)

#load data
d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

#choose which study
#dataSets <- unique(allYrs[,c("country","study_id")])
#write.table(dataSets,file="datasets.txt",sep="\t",row.names=FALSE)
TaskID <- read.delim("/data/idiv_ess/Ellen/datasets.txt",as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
mytask <- TaskID$study_id[which(TaskID$task_id==task.id)]
allYrs <- subset(allYrs,study_id==mytask)

#centre Year - helps model convergence to center variables for the model
allYrs$cYear <- allYrs$year_wMissing - median(allYrs$year_wMissing)
summary(allYrs$cYear)
#or just have as an index starting from 1
allYrs$iYear <- allYrs$year_wMissing - min(allYrs$year_wMissing)+1
summary(allYrs$iYear)

#centre day of year
allYrs$cday_of_year <- allYrs$day_of_year - median(allYrs$day_of_year,na.rm=T)

#### two-stage models ####

#### fit model using brms ####

# library(brms)
# #rstan_options(auto_write = TRUE)
# #options(mc.cores = parallel::detectCores())
# 
# #set priors now
# prior1 = c(set_prior("normal(0,1)", class = "b"))
# 
# # write a function to consider day of year if sampling is more than 30 days apart
# fitStanModel <- function(mydata){
#   
#   #if sampling occurs in more than one month include a seasonal term in the model
#   maxDiffDays = max(mydata$cday_of_year)-min(mydata$cday_of_year)
#   
#   if(maxDiffDays < 30) {
#     myformula <- bf(log10(spp_richness+1) ~ cYear + ar(time = iYear, p = 1, cov=FALSE))
#   } else{
#     myformula <- bf(log10(spp_richness+1) ~ cday_of_year + cYear + ar(time = iYear, p = 1, cov=FALSE))
#   }
#   
#   #fit model
#   fit1 <- brm(myformula, data = mydata, 
#               #family = poisson(), 
#               chains = n.cores,
#               prior = prior1, 
#               refresh = 0)
#   
#   #extract model fits
#   modelSummary <- fixef(fit1, pars="cYear")[1, c(1,2)]
#   return(modelSummary)
#   
# }
# 
# #get cores
# n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
# 
# #loop for all sites
# allsites <- sort(unique(allYrs$site_id))
# 
# trends <- lapply(allsites, function(x){
#   fitStanModel(subset(allYrs, site_id == x))
# })
# 
# trends <- data.frame(do.call(rbind, trends))
# trends$siteID <- allsites
# 
# saveRDS(trends, file=paste0("trends_",mytask,".RDS"))


### fitting directly in stan #####
library(brms)
library(rstan)

#set priors now
prior1 = c(set_prior("normal(0,10)", class = "b"))

#get model code
#model_code <- make_stancode(myformula, data = mydata, 
#                            #family = poisson(), 
#                            chains = n.cores,
#                            prior = prior1, 
#                            refresh = 0)

# write a function to consider day of year if sampling is more than 30 days apart
fitStanModel <- function(mydata){
  
  #if sampling occurs in more than one month include a seasonal term in the model
  maxDiffDays = max(mydata$cday_of_year)-min(mydata$cday_of_year)
  
  if(maxDiffDays < 30) {
    myformula <- bf(log10(spp_richness+1) ~ cYear + ar(time = iYear, p = 1, cov=FALSE))
    modelfile <- "/data/idiv_ess/Ellen/stan_code.stan"
    
  } else{
    myformula <- bf(log10(spp_richness+1) ~ cYear + cday_of_year + ar(time = iYear, p = 1, cov=FALSE))
    modelfile <- "/data/idiv_ess/Ellen/stan_code_seasonal.stan"
    
  }
  
  
  #get model data
  model_data <- make_standata(myformula, data = mydata, 
                              chains = n.cores,
                              refresh = 0)
  model_data$cYear <- mydata$cYear
  model_data$cday <- mydata$cday_of_year

  #fit model in stan
  stan_model <- stan(modelfile, 
                     data = model_data, seed = 20)
  
  
  #extract model fits
  modelSummary <- summary(stan_model)$summary
  modelFits <- data.frame(estimate = modelSummary[1,"mean"],
                             sd = modelSummary[1,"sd"])
  return(modelFits)
  
}

#get cores
n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

#loop for all sites
allsites <- sort(unique(allYrs$site_id))

trends <- lapply(allsites, function(x){
  fitStanModel(subset(allYrs, site_id == x))
})

trends <- data.frame(do.call(rbind, trends))
trends$siteID <- allsites

saveRDS(trends, file=paste0("trends_",mytask,".RDS"))
