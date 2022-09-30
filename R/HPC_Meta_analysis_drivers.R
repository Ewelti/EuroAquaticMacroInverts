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

#remove slopes for urban and crop and N

siteData <- unique(d1[,c("site","study_id","Country",
                         "ppt_Est", "tmax_Est", 
                         "ppt_mm_12moPrior", "tmax_C_12moPrior", 
                         "strahler_streamOrder","accumulation_atPoint", 
                         "elevation_atPoint","slope_mean", 
                         #"N_Est", "N_mean", 
                         "urban_meanPerc_upstr", "crop_meanPerc_upstr",
                         #"crop_Est", "urban_Est", 
                         "dam_impact_score_lessthan100km","TaxonomicRes")])

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
nrow(siteData)

### merge with biodiversity data ####

response_stan <- merge(siteData,response_stan,by="site_id")

### run model ####

library(brms)
library(rstan)

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

# #define priors 
# prior1 = c(set_prior("normal(0,1)", class = "b"))
# 
# #unweighted
# 
# fit1 <- brm(estimate ~ sppt_Est + stmax_Est + sppt_mm_12moPrior + stmax_C_12moPrior + 
#               sstrahler_streamOrder + saccumulation_atPoint + selevation_atPoint +
#               sslope_mean + surban_meanPerc_upstr + scrop_meanPerc_upstr +
#               sdam_impact_score_lessthan100km + (1|study_id) + (1|country),
#             data = response_stan, iter=5000, chains = 4, prior=prior1,
#             control = list(adapt_delta = 0.90, max_treedepth = 12))
# 
# #### save output ####
# 
# saveRDS(fit1,file=paste0("metaanalysis_unweighted_drivers_",myResponse,".rds"))
# 
# #weighted
# 
# fit1 <- brm(estimate|se(sd) ~ sppt_Est + stmax_Est + sppt_mm_12moPrior + stmax_C_12moPrior + 
#               sstrahler_streamOrder + saccumulation_atPoint + selevation_atPoint +
#               sslope_mean + surban_meanPerc_upstr + scrop_meanPerc_upstr +
#               sdam_impact_score_lessthan100km + (1|study_id) + (1|country),
#             data = response_stan, iter=5000, chains = 4, prior=prior1,
#             control = list(adapt_delta = 0.90, max_treedepth = 12))
# 
# #### save output ####
# 
# saveRDS(fit1,file=paste0("metaanalysis_drivers_",myResponse,".rds"))
# 
# ### prior check 1 ###
# 
# #also model with narrower prior - the horseshoe prior
# 
# prior1 = c(set_prior("horseshoe(1)", class = "b"))
# 
# fit1 <- brm(estimate|se(sd) ~ sppt_Est + stmax_Est + sppt_mm_12moPrior + stmax_C_12moPrior + 
#               sstrahler_streamOrder + saccumulation_atPoint + selevation_atPoint +
#               sslope_mean + surban_meanPerc_upstr + scrop_meanPerc_upstr +
#               sdam_impact_score_lessthan100km + (1|study_id) + (1|country),
#             data = response_stan, iter=5000, chains = 4, prior=prior1,
#             control = list(adapt_delta = 0.90, max_treedepth = 12))
# 
# #### save output ####
# 
# saveRDS(fit1,file=paste0("metaanalysis_drivers_horseshoe_",myResponse,".rds"))
# 
# ### prior check 2 ###
# 
# #also model with narrower prior - the lasso prior
# 
# prior1 = c(set_prior("lasso()", class = "b"))
# 
# fit1 <- brm(estimate|se(sd) ~ sppt_Est + stmax_Est + sppt_mm_12moPrior + stmax_C_12moPrior + 
#               sstrahler_streamOrder + saccumulation_atPoint + selevation_atPoint +
#               sslope_mean + surban_meanPerc_upstr + scrop_meanPerc_upstr +
#               sdam_impact_score_lessthan100km + (1|study_id) + (1|country),
#             data = response_stan, iter=5000, chains = 4, prior=prior1,
#             control = list(adapt_delta = 0.90, max_treedepth = 12))
# 
# #### save output ####
# 
# saveRDS(fit1,file=paste0("metaanalysis_drivers_lasso_",myResponse,".rds"))

### species only #####

prior1 = c(set_prior("horseshoe(1)", class = "b"))

fit1 <- brm(estimate|se(sd) ~ sppt_Est + stmax_Est + sppt_mm_12moPrior + stmax_C_12moPrior + 
              sstrahler_streamOrder + saccumulation_atPoint + selevation_atPoint +
              sslope_mean + surban_meanPerc_upstr + scrop_meanPerc_upstr +
              sdam_impact_score_lessthan100km + (1|study_id) + (1|country),
            data = subset(response_stan, TaxonomicRes=="species"), 
            iter=5000, chains = 4, prior=prior1,
            control = list(adapt_delta = 0.90, max_treedepth = 12))

saveRDS(fit1,file=paste0("metaanalysis_drivers_horseshoe_specieslevel_",myResponse,".rds"))
