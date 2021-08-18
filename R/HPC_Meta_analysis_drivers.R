#script to combine site-level trends together in a single meta-analysis and test the effect of the drivers

### get response for this task ######

TaskID <- read.csv("/data/idiv_ess/Ellen/ResponseTrends_TaskIDs.csv",as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]

### get site-level values for this response ####

response_stan <- readRDS("/data/idiv_ess/Ellen/stanTrends_site_level.rds")
response_stan <- subset(response_stan, Response == myResponse)
response_stan <- subset(response_stan, !is.na(estimate))

### site metadata ######

d1 <- read.csv("/data/idiv_ess/Ellen/All_siteLevel_and_glmOutput.csv", header=T) 

siteData <- unique(d1[,c("site","study_id","Country",
                         "ppt_Est", "tmax_Est", 
                         "ppt_mm_12moPrior", "tmax_C_12moPrior", 
                         "strahler_streamOrder","accumulation_atPoint", 
                         "elevation_atPoint","slope_mean", 
                         "N_Est", "N_mean", 
                         "urban_meanPerc_upstr", "crop_meanPerc_upstr",
                         "crop_Est", "urban_Est", 
                         "dam_impact_score_lessthan100km")])

#change names for consistency elsewhere
names(siteData)[which(names(siteData)=="site")] <- "site_id"
names(siteData)[which(names(siteData)=="Country")] <- "country"

### scale covariates ######

#function to add a new column onto the data with scaled vars (with s before their name)
scaleVars <- function(df){
  newd <- plyr::numcolwise(scale)(df)
  names(newd) <- sapply(names(newd),function(x)paste0("s",x))
  cbind(df, newd)
}

#apply function
siteData <- scaleVars(siteData)

### remove missing covariate data ####

siteData <- siteData[complete.cases(siteData),]
nrow(siteData)#1626

### merge with biodiversity data ####

response_stan <- merge(siteData,response_stan,by="site_id")

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

fit1 <- brm(estimate|weights(w) ~ sppt_Est + stmax_Est + sppt_mm_12moPrior + stmax_C_12moPrior + 
              sstrahler_streamOrder + saccumulation_atPoint + selevation_atPoint +
              sslope_mean + sN_Est + sN_mean + surban_meanPerc_upstr + scrop_meanPerc_upstr +
              scrop_Est + surban_Est + sdam_impact_score_lessthan100km + (1|study_id) + (1|country),
            data = response_stan, iter=4000, chains = n.chains)

#### save output ####

saveRDS(fit1,file=paste0("metaanalysis_drivers_",myResponse,".rds"))
