#script to combine site-level trends together in a single meta-analysis

### get response for this task ######

#we should restrict this to a few responses!!

TaskID <- read.csv("/data/idiv_ess/Ellen/ResponseTrends_TaskIDs.csv",as.is=T)
TaskID <- subset(TaskID,!duplicated(Response))
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]

### get site-level trends for this response ####

response_stan <- readRDS("/data/idiv_ess/Ellen/stanTrends_site_level.rds")
response_stan <- subset(response_stan, Response == myResponse)
response_stan <- subset(response_stan, !is.na(estimate))

### site metadata ######

d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T) 
d1<- d1[!is.na(d1$site_id_wMissing),]
siteData <- unique(d1[,c("site_id","study_id","country","season","TaxonomicRes")])
response_stan <- merge(siteData,response_stan,by="site_id")

### prepare for brms ####

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

#define priors - default mostly ok but set it for b
prior1 = c(set_prior("normal(0,10)", class = "b"))

### run season model ####

# #set spring as the reference level
# response_stan$season <- factor(response_stan$season,
#                                levels = c("spring","summer",
#                                           "fall","winter","No Match"))
# 
# #differences with respect to spring
# fit1 <- brm(estimate|weights(w) ~ season + (1|study_id) + (1|country),
#             data = response_stan, iter=4000, chains = 4,
#             prior = prior1)
# 
# #print output
# fixef(fit1)
# saveRDS(fixef(fit1), file=paste0("fixef_seasonDiff_",myResponse,".rds"))
# 
# #differences with respect to spring
# #unweighted
# fit1 <- brm(estimate ~ season + (1|study_id) + (1|country),
#             data = response_stan, iter=4000, chains = 4,
#             prior = prior1)
# 
# #print output
# fixef(fit1)
# saveRDS(fixef(fit1), file=paste0("fixef_seasonDiff_unweighted_",myResponse,".rds"))
# 
# #with intercept removed - each coefficient represents the mean trend per season
# fit1 <- brm(estimate|weights(w) ~ -1 + season + (1|study_id) + (1|country),
#             data = response_stan, iter=4000, chains = 4,
#             prior = prior1)          
#           
# saveRDS(fixef(fit1), file=paste0("fixef_seasonTrends_",myResponse,".rds"))
# 
# 
# #with intercept removed - each coefficient represents the mean trend per season
# #unweighted
# fit1 <- brm(estimate ~ -1 + season + (1|study_id) + (1|country),
#             data = response_stan, iter=4000, chains = 4,
#             prior = prior1)          
# 
# saveRDS(fixef(fit1), file=paste0("fixef_seasonTrends_unweighted_",myResponse,".rds"))

### run taxonomic resolution model ####

# #set species as the default level 
# response_stan$TaxonomicRes <- factor(response_stan$TaxonomicRes,
#                                      levels = c("species","genus","family"))
# 
# #difference in trends with respect to species
# fit1 <- brm(estimate|weights(w) ~ TaxonomicRes + (1|study_id) + (1|country),
#             data = response_stan, iter=4000, chains = 4,
#             prior = prior1)
# 
# #print output
# fixef(fit1)
# saveRDS(fixef(fit1), file=paste0("fixef_taxonresDiff_",myResponse,".rds"))
# 
# #difference in trends with respect to species
# #unweighted
# fit1 <- brm(estimate ~ TaxonomicRes + (1|study_id) + (1|country),
#             data = response_stan, iter=4000, chains = 4,
#             prior = prior1)
# 
# #print output
# fixef(fit1)
# saveRDS(fixef(fit1), file=paste0("fixef_taxonresDiff_unweighted",myResponse,".rds"))
# 
# 
# #trends for each taxonomic group
# fit1 <- brm(estimate|weights(w) ~ -1 + TaxonomicRes + (1|study_id) + (1|country),
#             data = response_stan, iter=4000, chains = 4,
#             prior = prior1)
# 
# #print output
# fixef(fit1)
# saveRDS(fixef(fit1), file=paste0("fixef_taxonresTrends_",myResponse,".rds"))
# 
# #trends for each taxonomic group
# #unweighted
# fit1 <- brm(estimate ~ -1 + TaxonomicRes + (1|study_id) + (1|country),
#             data = response_stan, iter=4000, chains = 4,
#             prior = prior1)
# 
# #print output
# fixef(fit1)
# saveRDS(fixef(fit1), file=paste0("fixef_taxonresTrends_unweighted_",myResponse,".rds"))

### taxonomic subsets #######

#weighted

#species
fit1_sp <- brm(estimate|weights(w) ~ 1 + (1|study_id) + (1|country),
            data = subset(response_stan, TaxonomicRes=="species"), 
            iter=4000, chains = 4,
            prior = prior1)
saveRDS(fixef(fit1_sp), file=paste0("fixef_taxonresSpecies_",myResponse,".rds"))

#genus
fit1_g <- brm(estimate|weights(w) ~ 1 + (1|study_id) + (1|country),
               data = subset(response_stan, TaxonomicRes=="genus"), 
               iter=4000, chains = 4,
               prior = prior1)
saveRDS(fixef(fit1_g), file=paste0("fixef_taxonresGenus_",myResponse,".rds"))

#family
fit1_f <- brm(estimate|weights(w) ~ 1 + (1|study_id) + (1|country),
              data = subset(response_stan, TaxonomicRes=="family"), 
              iter=4000, chains = 4,
              prior = prior1)
saveRDS(fixef(fit1_f), file=paste0("fixef_taxonresFamily_",myResponse,".rds"))

#as above but unweighted
#species
fit1_sp <- brm(estimate ~ 1 + (1|study_id) + (1|country),
               data = subset(response_stan, TaxonomicRes=="species"), 
               iter=4000, chains = 4,
               prior = prior1)
saveRDS(fixef(fit1_sp), file=paste0("fixef_taxonresSpecies_unweighted_",myResponse,".rds"))

#genus
fit1_g <- brm(estimate ~ 1 + (1|study_id) + (1|country),
              data = subset(response_stan, TaxonomicRes=="genus"), 
              iter=4000, chains = 4,
              prior = prior1)
saveRDS(fixef(fit1_g), file=paste0("fixef_taxonresGenus_unweighted_",myResponse,".rds"))

#family
fit1_f <- brm(estimate ~ 1 + (1|study_id) + (1|country),
              data = subset(response_stan, TaxonomicRes=="family"), 
              iter=4000, chains = 4,
              prior = prior1)
saveRDS(fixef(fit1_f), file=paste0("fixef_taxonresFamily_unweighted_",myResponse,".rds"))


