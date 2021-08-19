#read in model for each response and save fixed effects

### site-level trends ####
library(tidyverse)

response_stan <- readRDS("outputs/stanTrends_site_level.rds")
response_stan <- subset(response_stan, Response == "spp_richness")

#get site metadata
d1 <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears.csv", header=T)
d1<- d1[!is.na(d1$site_id_wMissing),]
siteData <- unique(d1[,c("site_id","study_id","country","season","TaxonomicRes")])
response_stan <- merge(siteData,response_stan,by="site_id")
head(response_stan)

summaryData <- response_stan %>%
  group_by(country,study_id,season,TaxonomicRes) %>%
  summarise(medTrends = median(estimate),
            nuData = length(estimate))

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
gls$stan_fit <- response_stan$estimate[match(gls$site,response_stan$site_id)]

qplot(stan_fit,SppRich_Est, data=gls)
cor(gls$stan_fit,gls$SppRich_Est)
#0.94

### meta-analysis ####

setwd("outputs/Meta-analysis")

library(brms)

#### spp_richness ####
fit <- readRDS("metaanalysis_spp_richness.rds")


#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
sr_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
sr_fixed <- list(Response="spp_richness", sr_fixed_99[,1:4], sr_fixed_975[,3:4],
                 sr_fixed_95[,3:4],sr_fixed_90[,3:4])
sr_fixed <-data.frame(lapply(sr_fixed, function(x) t(data.frame(x))))
sr_fixed 

#### spp_rich_rare ####
fit <- readRDS("metaanalysis_spp_rich_rare.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
srr_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
srr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
srr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
srr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
srr_fixed <- list(Response="spp_richness_rarefied", srr_fixed_99[,1:4], srr_fixed_975[,3:4],
                  srr_fixed_95[,3:4],srr_fixed_90[,3:4])
srr_fixed <-data.frame(lapply(srr_fixed, function(x) t(data.frame(x))))
srr_fixed

#### shannonsH ####
fit <- readRDS("metaanalysis_shannonsH.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
shH_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
shH_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
shH_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
shH_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
shH_fixed <- list(Response="shannonsH", shH_fixed_99[,1:4], shH_fixed_975[,3:4],
                  shH_fixed_95[,3:4],shH_fixed_90[,3:4])
shH_fixed <-data.frame(lapply(shH_fixed, function(x) t(data.frame(x))))
shH_fixed

#### E10 ####
fit <- readRDS("metaanalysis_E10.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
e10_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
e10_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
e10_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
e10_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
e10_fixed <- list(Response="E10", e10_fixed_99[,1:4], e10_fixed_975[,3:4],
                  e10_fixed_95[,3:4],e10_fixed_90[,3:4])
e10_fixed <-data.frame(lapply(e10_fixed, function(x) t(data.frame(x))))
e10_fixed

#### abundance ####
fit <- readRDS("metaanalysis_abundance.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
abund_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
abund_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
abund_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
abund_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
abund_fixed <- list(Response="abundance", abund_fixed_99[,1:4], abund_fixed_975[,3:4],
                    abund_fixed_95[,3:4],abund_fixed_90[,3:4])
abund_fixed <-data.frame(lapply(abund_fixed, function(x) t(data.frame(x))))
abund_fixed

#### turnover ####
fit <- readRDS("metaanalysis_turnover.rds") ## not yet calculated

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
turn_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
turn_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
turn_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
turn_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
turn_fixed <- list(Response="turnover", turn_fixed_99[,1:4], turn_fixed_975[,3:4],
                   turn_fixed_95[,3:4],turn_fixed_90[,3:4])
turn_fixed <-data.frame(lapply(turn_fixed, function(x) t(data.frame(x))))
turn_fixed

#### F_to ####
fit <- readRDS("metaanalysis_F_to.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
fto_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
fto_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fto_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fto_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fto_fixed <- list(Response="func_turnover", fto_fixed_99[,1:4], fto_fixed_975[,3:4],
                  fto_fixed_95[,3:4],fto_fixed_90[,3:4])
fto_fixed <-data.frame(lapply(fto_fixed, function(x) t(data.frame(x))))
fto_fixed

#### FRic ####
fit <- readRDS("metaanalysis_FRic.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
fric_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
fric_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fric_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fric_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fric_fixed <- list(Response="func_rich", fric_fixed_99[,1:4], fric_fixed_975[,3:4],
                   fric_fixed_95[,3:4],fric_fixed_90[,3:4])
fric_fixed <-data.frame(lapply(fric_fixed, function(x) t(data.frame(x))))
fric_fixed

#### FEve ####
fit <- readRDS("metaanalysis_FEve.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
feve_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
feve_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
feve_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
feve_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
feve_fixed <- list(Response="func_even", feve_fixed_99[,1:4], feve_fixed_975[,3:4],
                   feve_fixed_95[,3:4],feve_fixed_90[,3:4])
feve_fixed <-data.frame(lapply(feve_fixed, function(x) t(data.frame(x))))
feve_fixed

#### FDiv ####
fit <- readRDS("metaanalysis_FDiv.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
fdiv_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
fdiv_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fdiv_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fdiv_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fdiv_fixed <- list(Response="func_diverg", fdiv_fixed_99[,1:4], fdiv_fixed_975[,3:4],
                   fdiv_fixed_95[,3:4],fdiv_fixed_90[,3:4])
fdiv_fixed <-data.frame(lapply(fdiv_fixed, function(x) t(data.frame(x))))
fdiv_fixed

#### RaoQ ####
fit <- readRDS("metaanalysis_RaoQ.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
raoq_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
raoq_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
raoq_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
raoq_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
raoq_fixed <- list(Response="RaoQ", raoq_fixed_99[,1:4], raoq_fixed_975[,3:4],
                   raoq_fixed_95[,3:4],raoq_fixed_90[,3:4])
raoq_fixed <-data.frame(lapply(raoq_fixed, function(x) t(data.frame(x))))
raoq_fixed

#### alien_SppRich ####
fit <- readRDS("metaanalysis_alien_SppRich.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
aliensr_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
aliensr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
aliensr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
aliensr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
aliensr_fixed <- list(Response="alien_SppRich", aliensr_fixed_99[,1:4], aliensr_fixed_975[,3:4],
                      aliensr_fixed_95[,3:4],aliensr_fixed_90[,3:4])
aliensr_fixed <-data.frame(lapply(aliensr_fixed, function(x) t(data.frame(x))))
aliensr_fixed

#### alien_Abund ####
fit <- readRDS("metaanalysis_alien_Abund.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
alienab_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
alienab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
alienab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
alienab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
alienab_fixed <- list(Response="alien_abund", alienab_fixed_99[,1:4], alienab_fixed_975[,3:4],
                      alienab_fixed_95[,3:4],alienab_fixed_90[,3:4])
alienab_fixed <-data.frame(lapply(alienab_fixed, function(x) t(data.frame(x))))
alienab_fixed

#### abund_nativeSpp ####
fit <- readRDS("metaanalysis_abund_nativeSpp.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
nativeab_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
nativeab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
nativeab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
nativeab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
nativeab_fixed <- list(Response="native_abund", nativeab_fixed_99[,1:4], nativeab_fixed_975[,3:4],
                       nativeab_fixed_95[,3:4],nativeab_fixed_90[,3:4])
nativeab_fixed <-data.frame(lapply(nativeab_fixed, function(x) t(data.frame(x))))
nativeab_fixed

#### SppRich_nativeSpp ####
fit <- readRDS("metaanalysis_SppRich_nativeSpp.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
nativesr_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
nativesr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
nativesr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
nativesr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
nativesr_fixed <- list(Response="native_SppRich", nativesr_fixed_99[,1:4], nativesr_fixed_975[,3:4],
                       nativesr_fixed_95[,3:4],nativesr_fixed_90[,3:4])
nativesr_fixed <-data.frame(lapply(nativesr_fixed, function(x) t(data.frame(x))))
nativesr_fixed

#### EPT_SppRich ####
fit <- readRDS("metaanalysis_EPT_SppRich.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
EPTsr_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
EPTsr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
EPTsr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
EPTsr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
EPTsr_fixed <- list(Response="EPT_SppRich", EPTsr_fixed_99[,1:4], EPTsr_fixed_975[,3:4],
                    EPTsr_fixed_95[,3:4],EPTsr_fixed_90[,3:4])
EPTsr_fixed <-data.frame(lapply(EPTsr_fixed, function(x) t(data.frame(x))))
EPTsr_fixed

#### EPT_Abund ####
fit <- readRDS("metaanalysis_EPT_Abund.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
EPTab_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
EPTab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
EPTab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
EPTab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
EPTab_fixed <- list(Response="EPT_abund", EPTab_fixed_99[,1:4], EPTab_fixed_975[,3:4],
                    EPTab_fixed_95[,3:4],EPTab_fixed_90[,3:4])
EPTab_fixed <-data.frame(lapply(EPTab_fixed, function(x) t(data.frame(x))))
EPTab_fixed

#### insect_SppRich ####
fit <- readRDS("metaanalysis_insect_SppRich.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
insectsr_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
insectsr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
insectsr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
insectsr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
insectsr_fixed <- list(Response="insect_SppRich", insectsr_fixed_99[,1:4], insectsr_fixed_975[,3:4],
                       insectsr_fixed_95[,3:4],insectsr_fixed_90[,3:4])
insectsr_fixed <-data.frame(lapply(insectsr_fixed, function(x) t(data.frame(x))))
insectsr_fixed

#### insect_Abund ####
fit <- readRDS("metaanalysis_insect_Abund.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
insectab_fixed_99 <- fixef(fit, probs = c(0.01, 0.99))
insectab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
insectab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
insectab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
insectab_fixed <- list(Response="insect_abund", insectab_fixed_99[,1:4], insectab_fixed_975[,3:4],
                       insectab_fixed_95[,3:4],insectab_fixed_90[,3:4])
insectab_fixed <-data.frame(lapply(insectab_fixed, function(x) t(data.frame(x))))
insectab_fixed

#### assemble all model estimates from meta-analysis models #####

Yr_metaanaly_Ests <- rbind(sr_fixed, srr_fixed, shH_fixed, e10_fixed, abund_fixed, turn_fixed, 
                           fto_fixed, fric_fixed, feve_fixed, fdiv_fixed, raoq_fixed, aliensr_fixed,
                           alienab_fixed, nativesr_fixed, nativeab_fixed, EPTsr_fixed, EPTab_fixed,
                           insectsr_fixed, insectab_fixed)
write.csv(Yr_metaanaly_Ests, "Yr_metaanaly_Ests.csv")
