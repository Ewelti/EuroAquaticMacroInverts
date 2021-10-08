#script to combine site-level trends together in a single meta-analysis

### get response for this task ######
TaskID <- read.csv("/data/idiv_ess/Ellen/MovingAverage_Meta_TaskIDs.csv",as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

#and data
response_stan <- readRDS("/data/idiv_ess/Ellen/stanTrends_site_level_movingaverages.rds")

### get start year for this task #############

myStartYear <- TaskID$StartYear[which(TaskID$TaskID==task.id)]
response_stan <- subset(response_stan, StartYear == myStartYear)

### get site-level values for this response ####

myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]
response_stan <- subset(response_stan, Response == myResponse)
response_stan <- subset(response_stan, !is.na(estimate))
  
### site metadata ######

d1 <- read.csv("/data/idiv_ess/Ellen/All_siteLevel_and_glmOutput.csv", header=T) 
siteData <- unique(d1[,c("site","study_id","Country","Latitude_Y")])
names(siteData)[1] <- "site_id"
names(siteData)[3] <- "country"
medianLat <- median(siteData$Latitude_Y)
response_stan <- merge(siteData,response_stan,by="site_id")

### run model ####

library(rstan)
library(brms)

### decide on priors ####

prior1 = c(set_prior("normal(0,5)", class = "Intercept"))

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


#data above median Lat

#define priors
fit1 <- brm(estimate|weights(w) ~ 1 + (1|study_id) + (1|country),
            data = subset(response_stan,Latitude_Y>medianLat), 
            iter=3000, inits = 0,
            chains = 4, prior = prior1,
            init = "0",
            control = list(adapt_delta = 0.90, 
                           max_treedepth = 12))

### save output ####

saveRDS(fit1,file=paste0("metaanalysis_movingaverage_north_",myResponse,"_",myStartYear,".rds"))


#data below median Lat

#define priors
fit1 <- brm(estimate|weights(w) ~ 1 + (1|study_id) + (1|country),
            data = subset(response_stan,Latitude_Y<=medianLat), 
            iter=3000, inits = 0,
            chains = 4, prior = prior1,
            init = "0",
            control = list(adapt_delta = 0.90, 
                           max_treedepth = 12))

### save output ####

saveRDS(fit1,file=paste0("metaanalysis_movingaverage_south_",myResponse,"_",myStartYear,".rds"))