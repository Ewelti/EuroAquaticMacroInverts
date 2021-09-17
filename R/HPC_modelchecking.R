#### function to extract posterior distribution of the trends ####

getTrendProbability <- function(fit){
  mySamples <- posterior_samples(fit,pars="b_Intercept")
  data.frame(probIncrease = mean(mySamples>0),probDecrease = mean(mySamples<0))
}

### end of functions ############################################


### site-level trends ####

#read in model for each response and save fixed effects

library(tidyverse)

response_stan <- readRDS("outputs/stanTrends_site_level.rds")

#pivot responses
response_stan_pivot <- response_stan %>%
                    select(c(Response,estimate,site_id)) %>%
                    pivot_wider(names_from = "Response",
                                values_from = "estimate")

#get site metadata
d1 <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears.csv", header=T)
d1<- d1[!is.na(d1$site_id_wMissing),]
siteData <- unique(d1[,c("site_id","study_id","country","season","TaxonomicRes")])
response_stan_pivot <- merge(siteData,response_stan_pivot,by="site_id")
head(response_stan_pivot)

summaryData <- response_stan_pivot %>%
  group_by(country,study_id,season,TaxonomicRes) %>%
  summarise(medTrends = median(spp_richness),
            nuData = length(spp_richness))

ggplot(summaryData)+
  geom_text(aes(x=medTrends, y =nuData,label=study_id),size=2)+
  theme_classic()

ggplot(summaryData)+
  geom_boxplot(aes(x=TaxonomicRes, y =medTrends),size=2)+
  theme_classic()

ggplot(summaryData)+
  geom_boxplot(aes(x=season, y =medTrends),size=2)+
  theme_classic()


#check against gls fits
gls <- read.csv("outputs/All_siteLevel_and_glmOutput.csv",as.is=T)
head(all)
all <- merge(response_stan_pivot, gls, by.x="site_id", by.y="site")

#compare estimates
qplot(abund_nativeSpp, nativeAbun_Est, data=all)
qplot(abundance, Abun_Est, data=all)
qplot(alien_Abund, AlienAbun_Est, data=all)
qplot(alien_SppRich, AlienSppRich_Est, data=all)
qplot(insect_Abund, insect_Abund_Est, data=all)
qplot(insect_SppRich, insect_SppRich_Est, data=all)
qplot(EPT_Abund, EPT_Abund_Est, data=all)
qplot(EPT_SppRich, EPT_SppRich_Est, data=all)

qplot(F_to, F_to_Est, data=all)
qplot(FDiv, FDiv_Est, data=all)
qplot(FEve, FEve_Est, data=all)
qplot(FRic, FRic_Est, data=all)
qplot(RaoQ, RaoQ_Est, data=all)

qplot(shannonsH, ShanH_Est, data=all)
qplot(spp_rich_rare, SppRichRare_Est, data=all)
qplot(spp_richness, SppRich_Est, data=all)
qplot(SppRich_nativeSpp, nativeSppRich_Est, data=all)
qplot(turnover, TurnO_Est, data=all)

### meta-analysis ####

setwd("outputs/meta_updated")

#for unweighted models
setwd("C:/Users/db40fysa/Dropbox/Git/ellen_outputs/meta")


getwd()
library(rstan)
library(brms)
library(loo)

#### spp_richness ####
fit <- readRDS("metaanalysis_unweighted_spp_richness.rds")
fit <- readRDS("metaanalysis_spp_richness.rds")
loo_R2(fit)

#prob of trend
sr_prob <- getTrendProbability(fit)
sr_prob <- data.frame(Response="spp_richness", sr_prob[,1:2])
sr_prob

#check model
plot(fit)
sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
sr_loo
sr_parento <- as.list(pareto_k_table(sr_loo))
Count_sr <- rbind(sr_parento[[1]],sr_parento[[2]],sr_parento[[3]],sr_parento[[4]])
colnames(Count_sr) <- "SppRich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
sr_fixed <- list(Response="spp_richness", sr_fixed_995[,1:4], sr_fixed_975[,3:4],
                 sr_fixed_95[,3:4],sr_fixed_90[,3:4])
sr_fixed <-data.frame(lapply(sr_fixed, function(x) t(data.frame(x))))
sr_fixed 

#### spp_rich_rare ####
fit <- readRDS("metaanalysis_unweighted_spp_rich_rare.rds")
fit <- readRDS("metaanalysis_spp_rich_rare.rds")

#prob of trend
srr_prob <- getTrendProbability(fit)
srr_prob <- data.frame(Response="spp_richness_rarefied", srr_prob[,1:2])
srr_prob

#check model
plot(fit)
srr_loo <- loo(fit, cores = getOption("mc.cores", 1))
srr_loo
srr_parento <- as.list(pareto_k_table(srr_loo))
Count_srr <- rbind(srr_parento[[1]],srr_parento[[2]],srr_parento[[3]],srr_parento[[4]])
colnames(Count_srr) <- "SppRichRarefied"
pp_check(fit, nsamples = 100)

#pull out fixed effects
srr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
srr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
srr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
srr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
srr_fixed <- list(Response="spp_richness_rarefied", srr_fixed_995[,1:4], srr_fixed_975[,3:4],
                  srr_fixed_95[,3:4],srr_fixed_90[,3:4])
srr_fixed <-data.frame(lapply(srr_fixed, function(x) t(data.frame(x))))
srr_fixed

#### shannonsH ####
fit <- readRDS("metaanalysis_unweighted_shannonsH.rds")
fit <- readRDS("metaanalysis_shannonsH.rds")

#prob of trend
shH_prob <- getTrendProbability(fit)
shH_prob <- data.frame(Response="shannonsH", shH_prob[,1:2])
shH_prob

#check model
plot(fit)
shH_loo <- loo(fit, cores = getOption("mc.cores", 1))
shH_loo
shH_parento <- as.list(pareto_k_table(shH_loo))
Count_shH <- rbind(shH_parento[[1]],shH_parento[[2]],shH_parento[[3]],shH_parento[[4]])
colnames(Count_shH) <- "ShannonsH"
pp_check(fit, nsamples = 100)

#pull out fixed effects
shH_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
shH_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
shH_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
shH_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
shH_fixed <- list(Response="shannonsH", shH_fixed_995[,1:4], shH_fixed_975[,3:4],
                  shH_fixed_95[,3:4],shH_fixed_90[,3:4])
shH_fixed <-data.frame(lapply(shH_fixed, function(x) t(data.frame(x))))
shH_fixed

#### E10 ####
fit <- readRDS("metaanalysis_unweighted_E10.rds")
fit <- readRDS("metaanalysis_E10.rds")

#prob of trend
e10_prob <- getTrendProbability(fit)
e10_prob <- data.frame(Response="E10", e10_prob[,1:2])
e10_prob

#check model
plot(fit)
e10_loo <- loo(fit, cores = getOption("mc.cores", 1))
e10_loo
e10_parento <- as.list(pareto_k_table(e10_loo))
Count_e10 <- rbind(e10_parento[[1]],e10_parento[[2]],e10_parento[[3]],e10_parento[[4]])
colnames(Count_e10) <- "E10"
pp_check(fit, nsamples = 100)

#pull out fixed effects
e10_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
e10_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
e10_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
e10_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
e10_fixed <- list(Response="E10", e10_fixed_995[,1:4], e10_fixed_975[,3:4],
                  e10_fixed_95[,3:4],e10_fixed_90[,3:4])
e10_fixed <-data.frame(lapply(e10_fixed, function(x) t(data.frame(x))))
e10_fixed

#### abundance ####
fit <- readRDS("metaanalysis_unweighted_abundance.rds")
fit <- readRDS("metaanalysis_abundance.rds")

#prob of trend
ab_prob <- getTrendProbability(fit)
ab_prob <- data.frame(Response="abundance", ab_prob[,1:2])
ab_prob

#check model
plot(fit)
ab_loo <- loo(fit, cores = getOption("mc.cores", 1))
ab_loo
ab_parento <- as.list(pareto_k_table(ab_loo))
Count_ab <- rbind(ab_parento[[1]],ab_parento[[2]],ab_parento[[3]],ab_parento[[4]])
colnames(Count_ab) <- "Abundance"
pp_check(fit, nsamples = 100)

#pull out fixed effects
abund_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
abund_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
abund_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
abund_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
abund_fixed <- list(Response="abundance", abund_fixed_995[,1:4], abund_fixed_975[,3:4],
                    abund_fixed_95[,3:4],abund_fixed_90[,3:4])
abund_fixed <-data.frame(lapply(abund_fixed, function(x) t(data.frame(x))))
abund_fixed

#### turnover ####
fit <- readRDS("metaanalysis_unweighted_turnover.rds") 
fit <- readRDS("metaanalysis_turnover.rds") ## not yet calculated

#prob of trend
turn_prob <- getTrendProbability(fit)
turn_prob <- data.frame(Response="turnover", turn_prob[,1:2])
turn_prob

#check model
plot(fit)
turn_loo <- loo(fit, cores = getOption("mc.cores", 1))
turn_loo
turn_parento <- as.list(pareto_k_table(turn_loo))
Count_turn <- rbind(turn_parento[[1]],turn_parento[[2]],turn_parento[[3]],turn_parento[[4]])
colnames(Count_turn) <- "turnover"
pp_check(fit, nsamples = 100)

#pull out fixed effects
turn_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
turn_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
turn_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
turn_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
turn_fixed <- list(Response="turnover", turn_fixed_995[,1:4], turn_fixed_975[,3:4],
                   turn_fixed_95[,3:4],turn_fixed_90[,3:4])
turn_fixed <-data.frame(lapply(turn_fixed, function(x) t(data.frame(x))))
turn_fixed

#### F_to ####
fit <- readRDS("metaanalysis_unweighted_F_to.rds")
fit <- readRDS("metaanalysis_F_to.rds")

#prob of trend
fto_prob <- getTrendProbability(fit)
fto_prob <- data.frame(Response="func_turnover", fto_prob[,1:2])
fto_prob

#check model
plot(fit)#bad!!! - better when unweighted!!!
fto_loo <- loo(fit, cores = getOption("mc.cores", 1))
fto_loo
fto_parento <- as.list(pareto_k_table(fto_loo))
Count_fto <- rbind(fto_parento[[1]],fto_parento[[2]],fto_parento[[3]],fto_parento[[4]])
colnames(Count_fto) <- "func_turnover"
pp_check(fit, nsamples = 100)

#pull out fixed effects
fto_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
fto_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fto_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fto_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fto_fixed <- list(Response="func_turnover", fto_fixed_995[,1:4], fto_fixed_975[,3:4],
                  fto_fixed_95[,3:4],fto_fixed_90[,3:4])
fto_fixed <-data.frame(lapply(fto_fixed, function(x) t(data.frame(x))))
fto_fixed

#### FRic ####
fit <- readRDS("metaanalysis_unweighted_FRic.rds")
fit <- readRDS("metaanalysis_FRic.rds")

#prob of trend
fric_prob <- getTrendProbability(fit)
fric_prob <- data.frame(Response="func_rich", fric_prob[,1:2])
fric_prob

#check model
plot(fit)
fric_loo <- loo(fit, cores = getOption("mc.cores", 1))
fric_loo
fric_parento <- as.list(pareto_k_table(fric_loo))
Count_fric <- rbind(fric_parento[[1]],fric_parento[[2]],fric_parento[[3]],fric_parento[[4]])
colnames(Count_fric) <- "func_rich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
fric_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
fric_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fric_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fric_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fric_fixed <- list(Response="func_rich", fric_fixed_995[,1:4], fric_fixed_975[,3:4],
                   fric_fixed_95[,3:4],fric_fixed_90[,3:4])
fric_fixed <-data.frame(lapply(fric_fixed, function(x) t(data.frame(x))))
fric_fixed

#### FEve ####
fit <- readRDS("metaanalysis_unweighted_FEve.rds")
fit <- readRDS("metaanalysis_FEve.rds")

#prob of trend
feve_prob <- getTrendProbability(fit)
feve_prob <- data.frame(Response="func_even", feve_prob[,1:2])
feve_prob

#check model
plot(fit)
feve_loo <- loo(fit, cores = getOption("mc.cores", 1))
feve_loo
feve_parento <- as.list(pareto_k_table(feve_loo))
Count_feve <- rbind(feve_parento[[1]],feve_parento[[2]],feve_parento[[3]],feve_parento[[4]])
colnames(Count_feve) <- "func_even"
pp_check(fit, nsamples = 100)

#pull out fixed effects
feve_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
feve_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
feve_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
feve_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
feve_fixed <- list(Response="func_even", feve_fixed_995[,1:4], feve_fixed_975[,3:4],
                   feve_fixed_95[,3:4],feve_fixed_90[,3:4])
feve_fixed <-data.frame(lapply(feve_fixed, function(x) t(data.frame(x))))
feve_fixed

#### FDiv ####
fit <- readRDS("metaanalysis_unweighted_FDiv.rds")
fit <- readRDS("metaanalysis_FDiv.rds")

#prob of trend
fdiv_prob <- getTrendProbability(fit)
fdiv_prob <- data.frame(Response="func_diverg", fdiv_prob[,1:2])
fdiv_prob

#check model
plot(fit)
fdiv_loo <- loo(fit, cores = getOption("mc.cores", 1))
fdiv_loo
fdiv_parento <- as.list(pareto_k_table(fdiv_loo))
Count_fdiv <- rbind(fdiv_parento[[1]],fdiv_parento[[2]],fdiv_parento[[3]],fdiv_parento[[4]])
colnames(Count_fdiv) <- "func_diverg"
pp_check(fit, nsamples = 100)

#pull out fixed effects
fdiv_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
fdiv_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fdiv_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fdiv_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fdiv_fixed <- list(Response="func_diverg", fdiv_fixed_995[,1:4], fdiv_fixed_975[,3:4],
                   fdiv_fixed_95[,3:4],fdiv_fixed_90[,3:4])
fdiv_fixed <-data.frame(lapply(fdiv_fixed, function(x) t(data.frame(x))))
fdiv_fixed

#### RaoQ ####
fit <- readRDS("metaanalysis_unweighted_RaoQ.rds")
fit <- readRDS("metaanalysis_RaoQ.rds")

#prob of trend
raoq_prob <- getTrendProbability(fit)
raoq_prob <- data.frame(Response="RaoQ", raoq_prob[,1:2])
raoq_prob

#check model
plot(fit)
raoq_loo <- loo(fit, cores = getOption("mc.cores", 1))
raoq_loo
raoq_parento <- as.list(pareto_k_table(raoq_loo))
Count_raoq <- rbind(raoq_parento[[1]],raoq_parento[[2]],raoq_parento[[3]],raoq_parento[[4]])
colnames(Count_raoq) <- "RaoQ"
pp_check(fit, nsamples = 100)

#pull out fixed effects
raoq_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
raoq_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
raoq_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
raoq_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
raoq_fixed <- list(Response="RaoQ", raoq_fixed_995[,1:4], raoq_fixed_975[,3:4],
                   raoq_fixed_95[,3:4],raoq_fixed_90[,3:4])
raoq_fixed <-data.frame(lapply(raoq_fixed, function(x) t(data.frame(x))))
raoq_fixed

#### alien_SppRich ####
fit <- readRDS("metaanalysis_unweighted_alien_SppRich.rds")
fit <- readRDS("metaanalysis_alien_SppRich.rds")

#prob of trend
aliensr_prob <- getTrendProbability(fit)
aliensr_prob <- data.frame(Response="alien_sppRich", aliensr_prob[,1:2])
aliensr_prob

#check model
plot(fit)#bad!! - better unweighted
aliensr_loo <- loo(fit, cores = getOption("mc.cores", 1))
aliensr_loo
aliensr_parento <- as.list(pareto_k_table(aliensr_loo))
Count_aliensr <- rbind(aliensr_parento[[1]],aliensr_parento[[2]],aliensr_parento[[3]],aliensr_parento[[4]])
colnames(Count_aliensr) <- "alien_sppRich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
aliensr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
aliensr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
aliensr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
aliensr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
aliensr_fixed <- list(Response="alien_SppRich", aliensr_fixed_995[,1:4], aliensr_fixed_975[,3:4],
                      aliensr_fixed_95[,3:4],aliensr_fixed_90[,3:4])
aliensr_fixed <-data.frame(lapply(aliensr_fixed, function(x) t(data.frame(x))))
aliensr_fixed

#### alien_Abund ####
fit <- readRDS("metaanalysis_unweighted_alien_Abund.rds")
fit <- readRDS("metaanalysis_alien_Abund.rds")

#prob of trend
alienab_prob <- getTrendProbability(fit)
alienab_prob <- data.frame(Response="alien_abund", alienab_prob[,1:2])
alienab_prob

#check model
plot(fit) 
alienab_loo <- loo(fit, cores = getOption("mc.cores", 1))
alienab_loo
alienab_parento <- as.list(pareto_k_table(alienab_loo))
Count_alienab <- rbind(alienab_parento[[1]],alienab_parento[[2]],alienab_parento[[3]],alienab_parento[[4]])
colnames(Count_alienab) <- "alien_abund"
pp_check(fit, nsamples = 100)

#pull out fixed effects
alienab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
alienab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
alienab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
alienab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
alienab_fixed <- list(Response="alien_abund", alienab_fixed_995[,1:4], alienab_fixed_975[,3:4],
                      alienab_fixed_95[,3:4],alienab_fixed_90[,3:4])
alienab_fixed <-data.frame(lapply(alienab_fixed, function(x) t(data.frame(x))))
alienab_fixed

#### abund_nativeSpp ####
fit <- readRDS("metaanalysis_unweighted_abund_nativeSpp.rds")
fit <- readRDS("metaanalysis_abund_nativeSpp.rds")

#prob of trend
nativeab_prob <- getTrendProbability(fit)
nativeab_prob <- data.frame(Response="native_abund", nativeab_prob[,1:2])
nativeab_prob

#check model
plot(fit)# a bit weird - better unweighted
nativeab_loo <- loo(fit, cores = getOption("mc.cores", 1))
nativeab_loo
nativeab_parento <- as.list(pareto_k_table(nativeab_loo))
Count_nativeab <- rbind(nativeab_parento[[1]],nativeab_parento[[2]],nativeab_parento[[3]],nativeab_parento[[4]])
colnames(Count_nativeab) <- "native_abund"
pp_check(fit, nsamples = 100)

#pull out fixed effects
nativeab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
nativeab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
nativeab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
nativeab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
nativeab_fixed <- list(Response="native_abund", nativeab_fixed_995[,1:4], nativeab_fixed_975[,3:4],
                       nativeab_fixed_95[,3:4],nativeab_fixed_90[,3:4])
nativeab_fixed <-data.frame(lapply(nativeab_fixed, function(x) t(data.frame(x))))
nativeab_fixed

#### SppRich_nativeSpp ####
fit <- readRDS("metaanalysis_unweighted_SppRich_nativeSpp.rds")
fit <- readRDS("metaanalysis_SppRich_nativeSpp.rds")

#prob of trend
nativesr_prob <- getTrendProbability(fit)
nativesr_prob <- data.frame(Response="native_sppRich", nativesr_prob[,1:2])
nativesr_prob

#check model
plot(fit)
nativesr_loo <- loo(fit, cores = getOption("mc.cores", 1))
nativesr_loo
nativesr_parento <- as.list(pareto_k_table(nativesr_loo))
Count_nativesr <- rbind(nativesr_parento[[1]],nativesr_parento[[2]],nativesr_parento[[3]],nativesr_parento[[4]])
colnames(Count_nativesr) <- "native_sppRich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
nativesr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
nativesr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
nativesr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
nativesr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
nativesr_fixed <- list(Response="native_SppRich", nativesr_fixed_995[,1:4], nativesr_fixed_975[,3:4],
                       nativesr_fixed_95[,3:4],nativesr_fixed_90[,3:4])
nativesr_fixed <-data.frame(lapply(nativesr_fixed, function(x) t(data.frame(x))))
nativesr_fixed

#### EPT_SppRich ####
fit <- readRDS("metaanalysis_unweighted_EPT_SppRich.rds")
fit <- readRDS("metaanalysis_EPT_SppRich.rds")

#prob of trend
EPTsr_prob <- getTrendProbability(fit)
EPTsr_prob <- data.frame(Response="EPT_sppRich", EPTsr_prob[,1:2])
EPTsr_prob

#check model
plot(fit)#bad!! better unweighted
EPTsr_loo <- loo(fit, cores = getOption("mc.cores", 1))
EPTsr_loo
EPTsr_parento <- as.list(pareto_k_table(EPTsr_loo))
Count_EPTsr <- rbind(EPTsr_parento[[1]],EPTsr_parento[[2]],EPTsr_parento[[3]],EPTsr_parento[[4]])
colnames(Count_EPTsr) <- "EPT_sppRich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
EPTsr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
EPTsr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
EPTsr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
EPTsr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
EPTsr_fixed <- list(Response="EPT_SppRich", EPTsr_fixed_995[,1:4], EPTsr_fixed_975[,3:4],
                    EPTsr_fixed_95[,3:4],EPTsr_fixed_90[,3:4])
EPTsr_fixed <-data.frame(lapply(EPTsr_fixed, function(x) t(data.frame(x))))
EPTsr_fixed

#### EPT_Abund ####
fit <- readRDS("metaanalysis_unweighted_EPT_Abund.rds")
fit <- readRDS("metaanalysis_EPT_Abund.rds")

#prob of trend
EPTab_prob <- getTrendProbability(fit)
EPTab_prob <- data.frame(Response="EPT_abund", EPTab_prob[,1:2])
EPTab_prob

#check model
plot(fit)#bad!! better unweighted
EPTab_loo <- loo(fit, cores = getOption("mc.cores", 1))
EPTab_loo
EPTab_parento <- as.list(pareto_k_table(EPTab_loo))
Count_EPTab <- rbind(EPTab_parento[[1]],EPTab_parento[[2]],EPTab_parento[[3]],EPTab_parento[[4]])
colnames(Count_EPTab) <- "EPT_abund"
pp_check(fit, nsamples = 100)

#pull out fixed effects
EPTab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
EPTab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
EPTab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
EPTab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
EPTab_fixed <- list(Response="EPT_abund", EPTab_fixed_995[,1:4], EPTab_fixed_975[,3:4],
                    EPTab_fixed_95[,3:4],EPTab_fixed_90[,3:4])
EPTab_fixed <-data.frame(lapply(EPTab_fixed, function(x) t(data.frame(x))))
EPTab_fixed

#### insect_SppRich ####
fit <- readRDS("metaanalysis_unweighted_insect_SppRich.rds")
fit <- readRDS("metaanalysis_insect_SppRich.rds")

#prob of trend
insectsr_prob <- getTrendProbability(fit)
insectsr_prob <- data.frame(Response="insect_sppRich", insectsr_prob[,1:2])
insectsr_prob

#check model
plot(fit)
insectsr_loo <- loo(fit, cores = getOption("mc.cores", 1))
insectsr_loo
insectsr_parento <- as.list(pareto_k_table(insectsr_loo))
Count_insectsr <- rbind(insectsr_parento[[1]],insectsr_parento[[2]],insectsr_parento[[3]],insectsr_parento[[4]])
colnames(Count_insectsr) <- "insect_sppRich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
insectsr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
insectsr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
insectsr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
insectsr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
insectsr_fixed <- list(Response="insect_SppRich", insectsr_fixed_995[,1:4], insectsr_fixed_975[,3:4],
                       insectsr_fixed_95[,3:4],insectsr_fixed_90[,3:4])
insectsr_fixed <-data.frame(lapply(insectsr_fixed, function(x) t(data.frame(x))))
insectsr_fixed

#### insect_Abund ####
fit <- readRDS("metaanalysis_unweighted_insect_Abund.rds")
fit <- readRDS("metaanalysis_insect_Abund.rds")

#prob of trend
insectab_prob <- getTrendProbability(fit)
insectab_prob <- data.frame(Response="insect_abund", insectab_prob[,1:2])
insectab_prob

#check model
plot(fit)
insectab_loo <- loo(fit, cores = getOption("mc.cores", 1))
insectab_loo
insectab_parento <- as.list(pareto_k_table(insectab_loo))
Count_insectab <- rbind(insectab_parento[[1]],insectab_parento[[2]],insectab_parento[[3]],insectab_parento[[4]])
colnames(Count_insectab) <- "insect_abund"
pp_check(fit, nsamples = 100)

#pull out fixed effects
insectab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
insectab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
insectab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
insectab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
insectab_fixed <- list(Response="insect_abund", insectab_fixed_995[,1:4], insectab_fixed_975[,3:4],
                       insectab_fixed_95[,3:4],insectab_fixed_90[,3:4])
insectab_fixed <-data.frame(lapply(insectab_fixed, function(x) t(data.frame(x))))
insectab_fixed

#### assemble all model estimates from meta-analysis models #####

Yr_metaanaly_Ests <- rbind(sr_fixed, srr_fixed, shH_fixed, e10_fixed, abund_fixed, turn_fixed, 
                           fto_fixed, fric_fixed, feve_fixed, fdiv_fixed, raoq_fixed, aliensr_fixed,
                           alienab_fixed, nativesr_fixed, nativeab_fixed, EPTsr_fixed, EPTab_fixed,
                           insectsr_fixed, insectab_fixed)
write.csv(Yr_metaanaly_Ests, "Yr_metaanaly_Unweight_Ests.csv")

#### assemble all probabilities of increases/decreases from meta-analysis models #####

Yr_metaanaly_probs <- rbind(sr_prob, srr_prob, shH_prob, e10_prob, ab_prob, turn_prob, 
                            fto_prob, fric_prob, feve_prob, fdiv_prob, raoq_prob, aliensr_prob,
                            alienab_prob, nativesr_prob, nativeab_prob, EPTsr_prob, EPTab_prob,
                            insectsr_prob, insectab_prob)
write.csv(Yr_metaanaly_probs, "Yr_metaanaly_Unweight_probabilities.csv")

#### assemble model counts from Parento k diagnostic values from meta-analysis models #####

Yr_metaanaly_parento <- cbind(Count_sr, Count_srr, Count_shH, Count_e10, Count_ab, Count_turn, Count_fto,
                               Count_fric, Count_feve, Count_fdiv, Count_raoq, Count_aliensr,
                              Count_alienab, Count_nativesr, Count_nativeab, Count_EPTsr, Count_EPTab,
                              Count_insectsr, Count_insectab)
rownames(Yr_metaanaly_parento) <- c("good[-Inf, 0.5]","ok[0.5, 0.7]","bad[0.7, 1]","verybad[1, Inf]")
write.csv(Yr_metaanaly_parento, "Yr_meta_parento_Unweight_ModelCounts.csv")


