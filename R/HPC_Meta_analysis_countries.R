#script to combine site-level trends together in a single meta-analysis

### get id for this task ######

task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

### get site-level values ####

response_stan <- readRDS("/data/idiv_ess/Ellen/stanTrends_site_level.rds")

### site metadata ######

d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears.csv", header=T) 
siteData <- unique(d1[,c("site_id","study_id","country")])
response_stan <- merge(siteData,response_stan,by="site_id")

### choose country for this task

allCountries <- sort(unique(response_stan$country))#22
myCountry <- allCountries[task.id]
response_stan <- subset(response_stan, country == myCountry)

### run model ####

library(rstan)
library(brms)

### decide on priors ####

prior1 = c(set_prior("normal(0,1)", class = "Intercept"))

# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
rstan_options(auto_write = TRUE)
options(mc.cores = cpus_per_task)

#define weights
response_stan$w <- 1/response_stan$sd

### species richness ####

response_sta_sr <- subset(response_stan, Response=="spp_richness")

fit1 <- brm(estimate|weights(w) ~ 1 + (1|study_id),
            data = response_sta_sr, iter=5000, inits = 0,
            chains = 4, prior = prior1,
            init = "0",
            control = list(adapt_delta = 0.90, 
                           max_treedepth = 12))

### save output ####

saveRDS(fixef(fit1),file=paste0("metaanalysis__spp_richness__",myCountry,".rds"))

### abundance ####

response_sta_abund <- subset(response_stan, Response=="abundance")

fit1 <- brm(estimate|weights(w) ~ 1 + (1|study_id),
            data = response_sta_abund, iter=5000, inits = 0,
            chains = 4, prior = prior1,
            init = "0",
            control = list(adapt_delta = 0.90, 
                           max_treedepth = 12))

### save output ####

saveRDS(fixef(fit1),file=paste0("metaanalysis__abundance__",myCountry,".rds"))
