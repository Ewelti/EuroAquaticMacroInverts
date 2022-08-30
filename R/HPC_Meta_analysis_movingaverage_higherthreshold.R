#script to combine site-level trends together in a single meta-analysis

#get data
response_stan <- readRDS("/data/idiv_ess/Ellen/stanTrends_site_level_movingaveragesHTMV.rds")

### get response for this task ######
TaskID <- unique(response_stan[,c("StartYear","Response")])
TaskID$TaskID <- 1:nrow(TaskID)
nrow(TaskID)#135
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

### get start year for this task #############

myStartYear <- TaskID$StartYear[which(TaskID$TaskID==task.id)]
response_stan <- subset(response_stan, StartYear == myStartYear)

### get site-level values for this response ####

myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]
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

#define weights
response_stan$w <- 1/response_stan$sd
summary(response_stan$w)

# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
rstan_options(auto_write = FALSE)
options(mc.cores = cpus_per_task)

#define priors
fit1 <- brm(estimate|se(sd) ~ 1 + (1|study_id) + (1|country),
            data = response_stan, iter=3000, init = 0,
            chains = 4, prior = prior1,
            init = "0",
            control = list(adapt_delta = 0.90, 
                           max_treedepth = 12))

### save output ####
saveRDS(fit1,file=paste0("metaanalysis_movingaverage_higherthreshold",myResponse,"_",myStartYear,".rds"))