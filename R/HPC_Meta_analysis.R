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

d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears.csv", header=T) 
siteData <- unique(d1[,c("site_id","study_id","country")])
response_stan <- merge(siteData,response_stan,by="site_id")

### decide on priors ####


### run model ####

library(brms)

#examine response
#hist(response_stan$estimate)
summary(response_stan$estimate)

#define weights
response_stan$w <- 1/response_stan$sd
summary(response_stan$w)

n.chains = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

#define priors - default ok
#prior1 = c(set_prior("normal(0,10)", class = "Intercept"))

fit1 <- brm(estimate|weights(w) ~ 1 + (1|study_id) + (1|country),
            data = response_stan, iter=4000, inits = 0,
            chains = n.chains)

### save output ####

saveRDS(fit1,file=paste0("metaanalysis_",myResponse,".rds"))