#setwd("C:/Users/Ellen/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")
#
rm(list=ls())

library(lubridate)

##attach data
sites <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T) # change file name according to the time series to be analyzed
#delete missing data rows
DATA2 <- sites[!is.na(sites$site_id_wMissing),]

#make turnover numeric
DATA2$turnover <- as.numeric(DATA2$turnover)

#subset by year
sites_later <- subset(DATA2,year > 1999)
head(sites_later)

#count number of sampling years per site
yearcount <- aggregate(year ~ site_id, data = sites_later, FUN = length)
head(yearcount)

#subset for sites with more years
siteslater_long <- yearcount[which(yearcount$year > 14),]
head(siteslater_long)
nrow(siteslater_long)#556

#subset full dataset for sites with at least 15yrs of data in 2000 or later
allYrs <- subset(sites_later, site_id %in% siteslater_long$site_id)
head(allYrs)

#to create the TaskID's
timeWindow <- 10
minimumThreshold <- 6

SufficientSites <- lapply(2000:2014, function(x){
  allYrs2 <- subset(allYrs, year_wMissing >= x & year_wMissing < (x+timeWindow))
  siteSummary <- tapply(allYrs2$abundance,allYrs2$site_id,length)
  data.frame(StartYear = x, site_id = names(siteSummary)[siteSummary>=minimumThreshold])

})

SufficientSites <- do.call(rbind, SufficientSites)
SufficientSites$country <- allYrs$country[match(SufficientSites$site_id,allYrs$site_id)]
SufficientSites <- unique(SufficientSites[,c("StartYear","country")])
SufficientSites <- rbind(SufficientSites,SufficientSites,SufficientSites)
#SufficientSites$Response <- c(rep("E10",nrow(SufficientSites)/3),
#                              rep("FEve", nrow(SufficientSites)/3),
#                              rep("F_to", nrow(SufficientSites)/3))
#SufficientSites$Response <- c(rep("FRic",nrow(SufficientSites)/2),
#                             rep("FRed", nrow(SufficientSites)/2))
SufficientSites$Response <- c(rep("abundance",nrow(SufficientSites)/2),
                              rep("spp_richness", nrow(SufficientSites)/2))
SufficientSites$TaskID <- 1:nrow(SufficientSites)

#write.table(SufficientSites,"outputs/MovingAverageHigherThreshold_TaskIDs.csv",sep=",",row.names=FALSE)

#get task id
#TaskID <- read.csv("/data/idiv_ess/Ellen/MovingAverage_TaskIDs.csv",as.is=T)
TaskID <- SufficientSites
nrow(TaskID)#468
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
#"alien_SppRich","SppRich_nativeSpp","EPT_SppRich","insect_SppRich", Fred

#order by site site year
allYrs <- allYrs[order(allYrs$year_wMissing),]

#### subset by start year ####

#choose which start year for this task
StartYear <- TaskID$StartYear[which(TaskID$TaskID==task.id)]

#subset 10 year time span
allYrs <- subset(allYrs, year_wMissing>=StartYear & year_wMissing<(StartYear+timeWindow))

### select sites with enough data ###

allYrs <- subset(allYrs, !is.na(Response))
siteSummary <- tapply(allYrs$Response,allYrs$site_id,length)
sufficientSites <- names(siteSummary)[siteSummary>= minimumThreshold]
allYrs <- subset(allYrs, site_id %in% sufficientSites)                

nrow(allYrs)
### fitting directly in stan #####

library(brms)
library(rstan)

# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
rstan_options(auto_write = FALSE)
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

saveRDS(trends, file=paste0("trendsHTMV__",myResponse,"__",myCountry,"__",StartYear,".RDS"))







