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
siteData <- unique(d1[,c("site_id","study_id","country")])
response_stan <- merge(siteData,response_stan,by="site_id")

### run model ####

library(brms)

#examine response
hist(response_stan$estimate)
summary(response_stan$estimate)

#define weights
response_stan$w <- 1/response_stan$sd
summary(response_stan$w)

#function to add a new column onto the data with scaled vars (with s before their name)
scaleVars <- function(df){
  newd <- plyr::numcolwise(scale)(df)
  names(newd) <- sapply(names(newd),function(x)paste0("s",x))
  cbind(df, newd)
}
#apply function
response_stan <- scaleVars(response_stan)

#define priors - default ok
#prior1 = c(set_prior("normal(0,10)", class = "Intercept"))

n.chains = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

fit1 <- brm(estimate|weights(w) ~ sppt_Est + stmax_Est + sppt_mm_12moPrior + stmax_C_12moPrior + 
              sstrahler_streamOrder + saccumulation_atPoint + selevation_atPoint +
              sslope_mean + sN_Est + sN_mean + surban_meanPerc_upstr + scrop_meanPerc_upstr +
              scrop_Est + surban_Est + sdam_impact_score_lessthan100km +(1|study_id) + (1|country),
            data = response_stan, iter=4000, chains = n.chains)

#### save output ####

saveRDS(fit1,file=paste0("metaanalysis_",myResponse,".rds"))

#for later, pull out section for plotting
#fixed(fit1)
