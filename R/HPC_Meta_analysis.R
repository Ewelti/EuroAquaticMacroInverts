#script to combine site-level trends together in a single meta-analysis

### get response for this task ######

TaskID <- read.csv("/data/idiv_ess/Ellen/ResponseTrends_TaskIDs.csv",as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]

### get site-level values for this response ####

response_stan <- readRDS("/data/idiv_ess/Ellen/stanTrends_site_level.rds")
response_stan <- subset(response_stan, Response == myResponse)

### site metadata ######

d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears.csv", header=T) 
siteData <- unique(d1[,c("site_id","study_id","Country")])
response_stan <- merge(siteData,response_stan,by.x="site_id",by.y="siteID")

### run model ####

library(brms)

#examine response
hist(response_stan$estimate)
summary(response_stan$estimate)

#define weights
response_stan$w <- 1/response_stan$sd
summary(response_stan$w)

#define priors
prior1 = c(set_prior("normal(0,10)", class = "Intercept"))

fit1 <- brm(estimate|weights(w) ~ 1 + (1|study_id) + (1|Country),
            data = response_stan, prior = prior1, iter=4000)

#### save output ####

saveRDS(fit1,file=paste0("metaanalysis_",myResponse,".rds"))

#for later, pull out section for plotting
#fixed(fit1)
