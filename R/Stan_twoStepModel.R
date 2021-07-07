#install.packages("brms")
if (!require(devtools)) {
  install.packages("devtools")
}
library(brms)

#load data
d1 <- read.csv("raw-data/All_indices_benthicMacroInverts_AllYears.csv", header=T) # change file name according to the time series to be analyzed
allYrs <- d1[!is.na(d1$site_id_wMissing),]
head(allYrs)

#### two-stage models ####

#fit trend model to a single time-series
site100000001<-allYrs[which(allYrs$site_id=="100000001"),]
fit1 <- brm(spp_richness ~ year_wMissing, data = site100000001,family = poisson())
fit1

#see what the default priors are
get_prior(spp_richness ~ year_wMissing, data = site100000001, family = poisson())

#set priors now
prior = c(set_prior("normal(0,0.5)", class = "ar"), #I don't know what/if to change here
          set_prior("normal(0,5)", class = "b"))

#including autocorrelation (of the residuals)
fit1 <- brm(spp_richness ~ year_wMissing, autocor = cor_ar(~year_wMissing, p = 1),data = site100000001, family = poisson()) #'cor_brms' objects for 'autocor' is deprecated
fit1 <- brm(spp_richness ~ year_wMissing + ar(time = year_wMissing, p = 1),data = site100000001, family = poisson())

#cor_brms#including random effects
fit1 <- brm(spp_richness ~ year_wMissing + (1|year_wMissing) + ar(time = year_wMissing, p = 1),data = site100000001, family = poisson())

##load pre-calculated slopes
slopeData <- read.csv("outputs/All_siteLevel_and_glmOutput.csv", header=T) # change file name according to the time series to be analyzed
head(slopeData)

#combine trends in mixed model

prior1 <- prior(normal(0,5), class = b) +
  prior(cauchy(0,2), class = sd)

sr <- slopeData[!is.na(slopeData$SppRich_Est),]
fit1 <- brm(SppRich_Est|weights(SppRich_SE) ~ 1 + (1|site),data = sr, family = gaussian(), prior = prior1)
fit1 <- brm(SppRich_Est|se(SppRich_SE, sigma = TRUE) ~ 1 + (1|site),data = sr, family = gaussian(), prior = prior1) #Ellen made this and is unsure about it

#where w = 1/sd of the trend estimates

### one-stage model #####

prior1 = c(set_prior("normal(0,0.5)", class = "ar"),
           set_prior("normal(0,5)", class = "b"),
           set_prior("cauchy(0,2)", class = "sd"))


fit1 <- brm(Response ~ Year + (1 + Year|StudyID), 
            autocor = cor_ar(~Year, p = 1), 
            data = mydata, family = poisson())

#### model fits ####

ranef(fit1)

summary(fit1)

plot(fit1,ask = FALSE)

plot(conditional_effects(fit1),ask=FALSE)

#### model checking ####

loo(fit1)

pp_check(fit1)

WAIC(fit1)

#### model comparison ####

fit1 <- brm(rating ~ treat + period + carry,
            data = inhaler)
fit1 <- add_criterion(fit1, "waic")

fit2 <- brm(rating ~ treat + period + carry + (1|subject),
            data = inhaler)
fit2 <- add_criterion(fit2, "waic")

# compare both models
loo_compare(fit1, fit2, criterion = "waic")