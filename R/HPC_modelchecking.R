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
sr_fixed <- fixef(fit)

#### spp_rich_rare ####
fit <- readRDS("metaanalysis_spp_rich_rare.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
srr_fixed <-fixef(fit)

#### shannonsH ####
fit <- readRDS("metaanalysis_shannonsH.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
shH_fixed <- fixef(fit)

#### E10 ####
fit <- readRDS("metaanalysis_E10.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
e10_fixed <- fixef(fit)

#### abundance ####
fit <- readRDS("metaanalysis_abundance.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
abund_fixed <- fixef(fit)

#### turnover ####
fit <- readRDS("metaanalysis_turnover.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
turn_fixed <- fixef(fit)

#### F_to ####
fit <- readRDS("metaanalysis_F_to.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
fto_fixed <- fixef(fit)

#### FRic ####
fit <- readRDS("metaanalysis_FRic.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
fric_fixed <- fixef(fit)

#### FEve ####
fit <- readRDS("metaanalysis_FEve.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
feve_fixed <- fixef(fit)

#### FDiv ####
fit <- readRDS("metaanalysis_FDiv.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
fdiv <- fixef(fit)

#### RaoQ ####
fit <- readRDS("metaanalysis_RaoQ.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
raoq_fixed <- fixef(fit)

#### alien_SppRich ####
fit <- readRDS("metaanalysis_alien_SppRich.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
aliensr_fixed <- fixef(fit)

#### alien_Abund ####
fit <- readRDS("metaanalysis_alien_Abund.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
alienabund_fixed <- fixef(fit)

#### abund_nativeSpp ####
fit <- readRDS("metaanalysis_abund_nativeSpp.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
nativeabund_fixed <- fixef(fit)

#### SppRich_nativeSpp ####
fit <- readRDS("metaanalysis_SppRich_nativeSpp.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
nativesr_fixed <- fixef(fit)

#### EPT_SppRich ####
fit <- readRDS("metaanalysis_EPT_SppRich.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
eptsr_fixed <- fixef(fit)

#### EPT_Abund ####
fit <- readRDS("metaanalysis_EPT_Abund.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
eptabund_fixed <- fixef(fit)

#### insect_SppRich ####
fit <- readRDS("metaanalysis_insect_SppRich.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
insectsr_fixed <- fixef(fit)

#### insect_Abund ####
fit <- readRDS("metaanalysis_insect_Abund.rds")

#check model
plot(fit)
loo(fit, cores = getOption("mc.cores", 1))
pp_check(fit, nsamples = 100)

#pull out fixed effects
insectabund_fixed <- fixef(fit)
