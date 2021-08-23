#script to combine site-level trends together in a single meta-analysis

### get response for this task ######

#we should restrict this to a few responses!!

TaskID <- read.csv("/data/idiv_ess/Ellen/Sensitivity_TaskIDs.csv",as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
myResponse <- TaskID$Response[which(TaskID$TaskID==task.id)]

### get site-level trends for this response ####

response_stan <- readRDS("/data/idiv_ess/Ellen/stanTrends_site_level.rds")
response_stan <- subset(response_stan, Response == myResponse)
response_stan <- subset(response_stan, !is.na(estimate))

### site metadata ######

d1 <- read.csv("/data/idiv_ess/Ellen/All_indices_benthicMacroInverts_AllYears.csv", header=T) 
d1<- d1[!is.na(d1$site_id_wMissing),]
siteData <- unique(d1[,c("site_id","study_id","country","season","TaxonomicRes")])
response_stan <- merge(siteData,response_stan,by="site_id")

### prepare for brms ####

library(brms)

#examine response
#hist(response_stan$estimate)
summary(response_stan$estimate)

#define weights
response_stan$w <- 1/response_stan$sd
summary(response_stan$w)

n.chains = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

### run season model ####

#define priors - default mostly ok but set it for b
prior1 = c(set_prior("normal(0,10)", class = "b"))

#set spring as the reference level
response_stan$season <- factor(response_stan$season,
                               levels = c("spring","summer",
                                          "fall","winter","No Match"))

#differences with respect to spring
fit1 <- brm(estimate|weights(w) ~ season + (1|study_id) + (1|country),
            data = response_stan, iter=4000, chains = n.chains,
            prior = prior1)

#print output
fixef(fit1)
saveRDS(fixef(fit1), file=paste0("fixef_seasonDiff_",myResponse,".rds"))

#with intercept removed - each coefficient represents the mean trend per season
fit1 <- brm(estimate|weights(w) ~ -1 + season + (1|study_id) + (1|country),
            data = response_stan, iter=4000, chains = n.chains,
            prior = prior1)          
          
saveRDS(fixef(fit1), file=paste0("fixef_seasonTrends_",myResponse,".rds"))

### run taxonomic resolution model ####

#set species as the default level 
response_stan$TaxonomicRes <- factor(response_stan$TaxonomicRes,
                                     levels = c("species","genus","family"))

#difference in trends with respect to species
fit1 <- brm(estimate|weights(w) ~ TaxonomicRes + (1|study_id) + (1|country),
            data = response_stan, iter=4000, chains = n.chains,
            prior = prior1)

#print output
fixef(fit1)
saveRDS(fixef(fit1), file=paste0("fixef_taxonresDiff_",myResponse,".rds"))

#trends for each taxonomic group
fit1 <- brm(estimate|weights(w) ~ -1 + TaxonomicRes + (1|study_id) + (1|country),
            data = response_stan, iter=4000, chains = n.chains,
            prior = prior1)

#print output
fixef(fit1)
saveRDS(fixef(fit1), file=paste0("fixef_taxonresTrends_",myResponse,".rds"))

### spatial autocorrelation ####

#or use provided distance matrix??

# distance <- as.matrix(dist(response_stan[,c("Longitude_X","Latitude_Y")]))
# K <- nrow(sr)
# W <- array(0, c(K, K))
# W[distance == 1] <- 1 	
# 
# #lets fiddle with priors later
# 
# #compare models with and without spatial autocorrelation
# fit1 <- brm(estimate|weights(w) ~ 1 + (1|study_id) + car(W, type="icar"),
#             data = response_stan, data2 = list(W = W),
#             iter=4000, chains = n.chains)
# fit1 <- add_criterion(fit1, "waic")
# 
# fit2 <- brm(estimate|weights(w) ~ 1 + (1|study_id) + (1|country),
#             data = response_stan, data2 = list(W = W),
#             iter=4000, chains = n.chains)
# fit2 <- add_criterion(fit2, "waic")
# 
# loo(fit1, fit2, cores = getOption("mc.cores", 1))
# loo_compare(fit1, fit2, criterion = "waic")

### end ####