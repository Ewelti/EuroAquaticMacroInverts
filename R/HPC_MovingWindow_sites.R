rm(list=ls())

library(lubridate)

#load data
d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

#make turnover numeric
allYrs$turnover <- as.numeric(allYrs$turnover)

#first time to create the TaskID's
# SufficientSites <- lapply(1971:2011, function(x){
# 
#   timewindow <- 10
#   allYrs2 <- subset(allYrs, year_wMissing >= x & year_wMissing < (x+timewindow))
#   siteSummary <- tapply(allYrs2$abundance,allYrs2$site_id,length)
#   data.frame(StartYear = x, site_id = names(siteSummary)[siteSummary>=7])
# 
# })
# 
# SufficientSites <- do.call(rbind, SufficientSites)
# SufficientSites$country <- allYrs$country[match(SufficientSites$site_id,allYrs$site_id)]
# SufficientSites <- unique(SufficientSites[,c("StartYear","country")])
# SufficientSites <- rbind(SufficientSites,SufficientSites)
# SufficientSites$Response <- c(rep("abundance",nrow(SufficientSites)/2),
#                               rep("spp_richness", nrow(SufficientSites)/2))
# SufficientSites$TaskID <- 1:nrow(SufficientSites)
# write.table(SufficientSites,"outputs/MovingAverage_TaskIDs.csv",sep=",",row.names=FALSE)

#get task id
TaskID <- read.csv("/data/idiv_ess/Ellen/MovingAverage_TaskIDs2.csv",as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

### country ###
myCountry <- TaskID$country[which(TaskID$TaskID==task.id)]
allYrs <- subset(allYrs,country==myCountry)

### response ###

#choose which response for this task
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]
allYrs$Response <- allYrs[,myResponse]

#third run transformation
if(myResponse %in% c("abundance","alien_Abund","abund_nativeSpp",
                     "EPT_Abund","insect_Abund","FRic")){
  
  allYrs$Response <- log10(allYrs$Response+1) 
  
}else if(myResponse %in% c("E10","F_to")){
  
  allYrs$Response <- log10(allYrs$Response+0.01) 
  
}else if(myResponse %in% c("FDiv")){
  
  allYrs$Response <- allYrs$Response^2
  
}

#not transformed: shannonH, RaoQ, turnover, FEve, "spp_richness","spp_rich_rare",
#"alien_SppRich","SppRich_nativeSpp","EPT_SppRich","insect_SppRich"

#order by site site year
allYrs <- allYrs[order(allYrs$year_wMissing),]

#### subset by start year ####

#choose which start year for this task
StartYear <- TaskID$StartYear[which(TaskID$TaskID==task.id)]

#subset 10 year time span
timewindow <- 10
allYrs <- subset(allYrs, year_wMissing>=StartYear & year_wMissing<(StartYear+timewindow))

### select sites with enough data ###

allYrs <- subset(allYrs, !is.na(Response))

siteSummary <- tapply(allYrs$Response,allYrs$site_id,length)
sufficientSites <- names(siteSummary)[siteSummary>=7]
allYrs <- subset(allYrs, site_id %in% sufficientSites)                

### fitting directly in stan #####

library(brms)
library(rstan)

# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
#rstan_options(auto_write = TRUE)
options(mc.cores = cpus_per_task)

# write a function to consider day of year if sampling is more than 30 days apart
fitStanModel <- function(mydata){
  
  #centre Year - helps model convergence to center variables for the model
  mydata$cYear <- mydata$year_wMissing - median(mydata$year_wMissing)
  mydata$iYear <- mydata$year_wMissing - min(mydata$year_wMissing)+1
  
  #scale day of year
  mydata$cday_of_year <- (mydata$day_of_year - mean(mydata$day_of_year))/sd(mydata$day_of_year)
  
  #if sampling occurs in more than one month include a seasonal term in the model
  maxDiffDays = max(mydata$day_of_year)-min(mydata$day_of_year)
  
  if(maxDiffDays < 30) {
    myformula <- bf(Response ~ cYear + ar(time = iYear, p = 1, cov=FALSE))
    modelfile <- "/data/idiv_ess/Ellen/stan_code.stan"
    
  } else{
    myformula <- bf(Response ~ cYear + cday_of_year + ar(time = iYear, p = 1, cov=FALSE))
    modelfile <- "/data/idiv_ess/Ellen/stan_code_seasonal.stan"
    
  }
  
  #there is missing data, don't run the model
  if(all(is.na(mydata$Response))){
    
    modelFits <- data.frame(estimate = NA,
                            sd = NA,
                            rhat = NA)
  }else{
    
  #get model data
  
  #remove any missing values
  mydata <- subset(mydata, !is.na(Response))
      
  model_data <- make_standata(myformula, data = mydata, 
                                chains = 4)
  model_data$cYear <- mydata$cYear
  model_data$cday <- mydata$cday_of_year
  model_data$meanResponse <- round(median(mydata$Response), 1)
  model_data$sdResponse <- max(round(mad(mydata$Response), 1), 2.5)
  
  #fit model in stan
  stan_model <- stan(modelfile, 
                     data = model_data, 
                     chains = 4,
                     iter = 3000)
  
  #extract model fits
  modelSummary <- summary(stan_model)$summary
  modelFits <- data.frame(estimate = modelSummary[1,"mean"],
                          sd = modelSummary[1,"sd"],
                          rhat = modelSummary[1,"Rhat"])
  
  }
  return(modelFits)
  
}


#loop for all sites
allsites <- sort(unique(allYrs$site_id))

trends <- lapply(allsites, function(x){
  fitStanModel(subset(allYrs, site_id == x))
})

trends <- data.frame(do.call(rbind, trends))
trends$siteID <- allsites

saveRDS(trends, file=paste0("trends__",myResponse,"__",myCountry,"__",StartYear,".RDS"))
