rm(list=ls())

library(lubridate)

#load data
d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T) 
allYrs <- d1[!is.na(d1$site_id_wMissing),]

#make turnover numeric
allYrs$turnover <- as.numeric(allYrs$turnover)

#choose which country for this task
TaskID <- read.csv("/data/idiv_ess/Ellen/ResponseTrends_TaskIDs.csv",as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

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

#still problems with func_turn(F_to), alien_SppRich, EPT_SppRich, and EPT_abund

#order by site site year
allYrs <- allYrs[order(allYrs$year_wMissing),]

#### two-stage models ####

### fitting directly in stan #####

library(brms)
library(rstan)
 
mydata <- allYrs

#centre Year - helps model convergence to center variables for the model
mydata$cYear <- mydata$year_wMissing - median(mydata$year_wMissing)

#or just have as an index starting from 1
mydata$iYear <- mydata$year_wMissing - min(mydata$year_wMissing)+1
  
#scale day of year
mydata$cday_of_year <- (mydata$day_of_year - mean(mydata$day_of_year))/sd(mydata$day_of_year)
  
# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
rstan_options(auto_write = TRUE)
options(mc.cores = cpus_per_task)

#decide on priors
prior1 = c(set_prior("normal(0,5)", class = "Intercept"))

fit1 <- brm(Response ~ cYear + cday_of_year + 
              (1 + cYear|site_id) + 
              (1 + cYear|study_id) + 
              (1 + cYear|country) + 
              ar(time = iYear, gr = site_id, p = 1, cov=FALSE),
            data = mydata, 
            iter=5000, 
            thin=2,
            chains = 4, 
            prior = prior1,
            control = list(adapt_delta = 0.90, max_treedepth = 12))

### save output ####
saveRDS(fit1,file=paste0("onestage_",myResponse,".rds"))