#script to combine site-level trends together in a single meta-analysis

### get response for this task ######

TaskID <- read.csv("/data/idiv_ess/Ellen/ResponseTrends_TaskIDs.csv",as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]

### get site-level values for this response ####

response_stan <- readRDS("/data/idiv_ess/Ellen/stanTrends_site_level.rds")
response_stan <- subset(response_stan, Response == myResponse)
response_stan <- subset(response_stan, !is.na(estimate))
  
### site metadata ######

d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T) 
siteData <- unique(d1[,c("site_id","study_id","country")])
response_stan <- merge(siteData,response_stan,by="site_id")

### run model ####

library(rstan)
library(brms)

### decide on priors ####

prior1 = c(set_prior("normal(0,3)", class = "Intercept"))

#examine response
#hist(response_stan$estimate)
summary(response_stan$estimate)

# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
rstan_options(auto_write = FALSE)
options(mc.cores = cpus_per_task)

#fit weighted model - removing data for each country one by one

allCountries <- sort(unique(siteData$country))

modelOutputs <- data.frame()

for(i in 1:length(allCountries)){

fit1 <- brm(estimate|se(sd) ~ 1 + (1|study_id) + (1|country),
            data = subset(response_stan, country!=allCountries[i]), 
            iter=5000, init = 0,
            chains = 4, prior = prior1)

  out <- summary(fit1)$fixed
  out$country <- allCountries[i]
  
  modelOutputs <- rbind(modelOutputs,out)
}

### save output ####
saveRDS(modelOutputs,file=paste0("metaanalysis_country_jacknife_",myResponse,".rds"))
summary(modelOutputs)
