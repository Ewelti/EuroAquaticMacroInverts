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
  newd <- lapply(df, function(x) if(is.numeric(x)){
    scale(x, center=TRUE, scale=TRUE)
  } else x)
  names(newd) <- sapply(names(newd),function(x)paste0("s",x))
  cbind(df[,c(1:3,15)], newd)
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
### prior check 1 ###

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

# ### species only #####
# 
# prior1 = c(set_prior("horseshoe(1)", class = "b"))
# 
# fit1 <- brm(estimate|se(sd) ~ sppt_Est + stmax_Est + sppt_mm_12moPrior + stmax_C_12moPrior + 
#               sstrahler_streamOrder + saccumulation_atPoint + selevation_atPoint +
#               sslope_mean + surban_meanPerc_upstr + scrop_meanPerc_upstr +
#               sdam_impact_score_lessthan100km + (1|study_id) + (1|country),
#             data = subset(response_stan, TaxonomicRes=="species"), 
#             iter=5000, chains = 4, prior=prior1,
#             control = list(adapt_delta = 0.90, max_treedepth = 12))
# 
# saveRDS(fit1,file=paste0("metaanalysis_drivers_horseshoe_specieslevel_",myResponse,".rds"))

### test ####################

prior1 = c(set_prior("horseshoe(1)", class = "b"))

fit1 <- brm(estimate ~ sppt_Est + stmax_Est + sppt_mm_12moPrior + stmax_C_12moPrior +
              (1|study_id) + (1|country),
            data = response_stan, iter=1000, prior=prior1)

#### All 11 driver variables: sppt_Est + stmax_Est + sppt_mm_12moPrior + stmax_C_12moPrior + 
#               sstrahler_streamOrder + saccumulation_atPoint + selevation_atPoint +
#               sslope_mean + surban_meanPerc_upstr + scrop_meanPerc_upstr +
#               sdam_impact_score_lessthan100km

#make table for prediction for precipitation trend effect, controlling for others at medians
library(marginaleffects)
df <- datagrid(newdata = response_stan, 
               FUN_numeric = median,
               sppt_Est = sort(unique(response_stan$sppt_Est)))
head(df)

#get prediction of the model for this set of covariate values
preds <- posterior_linpred(fit1, newdata = df, re_formula = NA) 
df$trend <- apply(preds, 2, mean) 
df$lower <- apply(preds, 2, function(x) quantile(x,0.025)) 
df$upper <- apply(preds, 2, function(x) quantile(x,0.975)) 

saveRDS(df, file="preds_sppt_Est.rds")

#make table for prediction for temperature trend effect, controlling for others at medians
library(marginaleffects)
df <- datagrid(newdata = response_stan, 
               FUN_numeric = median,
               stmax_Est = sort(unique(response_stan$stmax_Est)))
head(df)

#get prediction of the model for this set of covariate values
preds <- posterior_linpred(fit1, newdata = df, re_formula = NA) 
df$trend <- apply(preds, 2, mean) 
df$lower <- apply(preds, 2, function(x) quantile(x,0.025)) 
df$upper <- apply(preds, 2, function(x) quantile(x,0.975)) 

saveRDS(df, file="preds_stmax_Est.rds")

#make table for prediction for precipitation mean effect, controlling for others at medians
library(marginaleffects)
df <- datagrid(newdata = response_stan, 
               FUN_numeric = median,
               sppt_mm_12moPrior = sort(unique(response_stan$sppt_mm_12moPrior)))
head(df)

#get prediction of the model for this set of covariate values
preds <- posterior_linpred(fit1, newdata = df, re_formula = NA) 
df$trend <- apply(preds, 2, mean) 
df$lower <- apply(preds, 2, function(x) quantile(x,0.025)) 
df$upper <- apply(preds, 2, function(x) quantile(x,0.975)) 

saveRDS(df, file="preds_sppt_mm_12moPrior.rds")

#make table for prediction for temperature mean effect, controlling for others at medians
library(marginaleffects)
df <- datagrid(newdata = response_stan, 
               FUN_numeric = median,
               stmax_C_12moPrior = sort(unique(response_stan$stmax_C_12moPrior)))
head(df)

#get prediction of the model for this set of covariate values
preds <- posterior_linpred(fit1, newdata = df, re_formula = NA) 
df$trend <- apply(preds, 2, mean) 
df$lower <- apply(preds, 2, function(x) quantile(x,0.025)) 
df$upper <- apply(preds, 2, function(x) quantile(x,0.975)) 

saveRDS(df, file="preds_stmax_C_12moPrior.rds")

#make table for prediction for stream order effect, controlling for others at medians
library(marginaleffects)
df <- datagrid(newdata = response_stan, 
               FUN_numeric = median,
               sstrahler_streamOrder = sort(unique(response_stan$sstrahler_streamOrder)))
head(df)

#get prediction of the model for this set of covariate values
preds <- posterior_linpred(fit1, newdata = df, re_formula = NA) 
df$trend <- apply(preds, 2, mean) 
df$lower <- apply(preds, 2, function(x) quantile(x,0.025)) 
df$upper <- apply(preds, 2, function(x) quantile(x,0.975)) 

saveRDS(df, file="preds_sstrahler_streamOrder.rds")

#make table for prediction for accumulation effect, controlling for others at medians
library(marginaleffects)
df <- datagrid(newdata = response_stan, 
               FUN_numeric = median,
               saccumulation_atPoint = sort(unique(response_stan$saccumulation_atPoint)))
head(df)

#get prediction of the model for this set of covariate values
preds <- posterior_linpred(fit1, newdata = df, re_formula = NA) 
df$trend <- apply(preds, 2, mean) 
df$lower <- apply(preds, 2, function(x) quantile(x,0.025)) 
df$upper <- apply(preds, 2, function(x) quantile(x,0.975)) 

saveRDS(df, file="preds_saccumulation_atPoint.rds")

#make table for prediction for elevation effect, controlling for others at medians
library(marginaleffects)
df <- datagrid(newdata = response_stan, 
               FUN_numeric = median,
               selevation_atPoint = sort(unique(response_stan$selevation_atPoint)))
head(df)

#get prediction of the model for this set of covariate values
preds <- posterior_linpred(fit1, newdata = df, re_formula = NA) 
df$trend <- apply(preds, 2, mean) 
df$lower <- apply(preds, 2, function(x) quantile(x,0.025)) 
df$upper <- apply(preds, 2, function(x) quantile(x,0.975)) 

saveRDS(df, file="preds_selevation_atPoint.rds")

#make table for prediction for slope effect, controlling for others at medians
library(marginaleffects)
df <- datagrid(newdata = response_stan, 
               FUN_numeric = median,
               sslope_mean = sort(unique(response_stan$sslope_mean)))
head(df)

#get prediction of the model for this set of covariate values
preds <- posterior_linpred(fit1, newdata = df, re_formula = NA) 
df$trend <- apply(preds, 2, mean) 
df$lower <- apply(preds, 2, function(x) quantile(x,0.025)) 
df$upper <- apply(preds, 2, function(x) quantile(x,0.975)) 

saveRDS(df, file="preds_sslope_mean.rds")

#make table for prediction for urban cover effect, controlling for others at medians
library(marginaleffects)
df <- datagrid(newdata = response_stan, 
               FUN_numeric = median,
               surban_meanPerc_upstr = sort(unique(response_stan$surban_meanPerc_upstr)))
head(df)

#get prediction of the model for this set of covariate values
preds <- posterior_linpred(fit1, newdata = df, re_formula = NA) 
df$trend <- apply(preds, 2, mean) 
df$lower <- apply(preds, 2, function(x) quantile(x,0.025)) 
df$upper <- apply(preds, 2, function(x) quantile(x,0.975)) 

saveRDS(df, file="preds_surban_meanPerc_upstr.rds")

#make table for prediction for crop cover effect, controlling for others at medians
library(marginaleffects)
df <- datagrid(newdata = response_stan, 
               FUN_numeric = median,
               scrop_meanPerc_upstr = sort(unique(response_stan$scrop_meanPerc_upstr)))
head(df)

#get prediction of the model for this set of covariate values
preds <- posterior_linpred(fit1, newdata = df, re_formula = NA) 
df$trend <- apply(preds, 2, mean) 
df$lower <- apply(preds, 2, function(x) quantile(x,0.025)) 
df$upper <- apply(preds, 2, function(x) quantile(x,0.975)) 

saveRDS(df, file="preds_scrop_meanPerc_upstr.rds")
              
#make table for prediction for precipitation trend effect, controlling for others at medians
library(marginaleffects)
df <- datagrid(newdata = response_stan, 
               FUN_numeric = median,
               sdam_impact_score_lessthan100km = sort(unique(response_stan$sdam_impact_score_lessthan100km)))
head(df)

#get prediction of the model for this set of covariate values
preds <- posterior_linpred(fit1, newdata = df, re_formula = NA) 
df$trend <- apply(preds, 2, mean) 
df$lower <- apply(preds, 2, function(x) quantile(x,0.025)) 
df$upper <- apply(preds, 2, function(x) quantile(x,0.975)) 

saveRDS(df, file="preds_sdam_impact_score_lessthan100km.rds")


#quick plot (we wont run this on HPC)
ggplot(df) +
  geom_line(aes(x=sppt_Est, y=trend)) +
  geom_ribbon(aes(x=sppt_Est, ymin=lower, ymax=upper),alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed")

### end of test
