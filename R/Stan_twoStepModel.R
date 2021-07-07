#install.packages("brms")
# if (!require(devtools)) {
#   install.packages("devtools")
# }
library(brms)
library(lubridate)

#load data
d1 <- read.csv("raw-data/All_indices_benthicMacroInverts_AllYears.csv", header=T) # change file name according to the time series to be analyzed
allYrs <- d1[!is.na(d1$site_id_wMissing),]
head(allYrs)

#centre Year - helps model convergence to center variables for the model
allYrs$cYear <- allYrs$year_wMissing - median(allYrs$year_wMissing)
summary(allYrs$cYear)
#or just have as an index starting from 1
allYrs$iYear <- allYrs$year_wMissing - min(allYrs$year_wMissing)+1
summary(allYrs$iYear)

#pull out date information (check with Ellen that is right?)
allYrs$Date <- sapply(as.character(allYrs$sample_id),function(x)strsplit(x,"_")[[1]][2])
allYrs$Date <- as.Date(allYrs$Date, format="%d.%m.%Y")
allYrs$day_of_year <- yday(allYrs$Date)
allYrs$cday_of_year <- allYrs$day_of_year - median(allYrs$day_of_year,na.rm=T)
  
#### two-stage models ####

#fit trend model to a single time-series
site100000001<-allYrs[which(allYrs$site_id=="100000001"),]

#simplest model:
fit1 <- brm(spp_richness ~ year_wMissing, data = site100000001,family = poisson())
fit1

#see what the default priors are
get_prior(spp_richness ~ year_wMissing, data = site100000001, family = poisson())

#set priors now
prior1 = c(set_prior("normal(0,0.5)", class = "ar"), #I don't know what/if to change here
          set_prior("normal(0,5)", class = "b"))

#including autocorrelation (of the residuals)
#fit1 <- brm(spp_richness ~ year_wMissing, autocor = cor_ar(~year_wMissing, p = 1),data = site100000001, family = poisson()) #'cor_brms' objects for 'autocor' is deprecated

fit1 <- brm(spp_richness ~ year_wMissing + ar(time = year_wMissing, p = 1),data = site100000001, family = poisson())
summary(fit1)

#including seasonal term??
unique(site100000001$month)

fit1 <- brm(spp_richness ~ cday_of_year + cYear + ar(time = iYear, p = 1),data = site100000001, family = poisson(), prior = prior1)
summary(fit1)

#might need to consider month (or day of year) as a spline term though for other datasets??
#the above linear term would not model complex seasonal patterns...
fit1 <- brm(spp_richness ~ s(cday_of_year) + cYear + ar(time = iYear, p = 1),data = site100000001, family = poisson(), prior = prior1)

# write a function to do this for any given dataset
fitStanModel <- function(mydata){
  
  #write model formula
  
  #if sampling occurs in more than one month include a seasonal term in the model
  nuMonths = length(unique(mydata$month))
  
  if(nuMonths == 1) {
    myformula <- bf(spp_richness ~ cYear + ar(time = iYear, p = 1))
  } else{
    myformula <- bf(spp_richness ~ cday_of_year + cYear + ar(time = iYear, p = 1))
  }
  
  #fit model
  fit1 <- brm(myformula, data = mydata, family = poisson(), prior = prior1, refresh = 0)
  
  #extract model fits
  modelSummary <- fixef(fit1, pars="cYear")
  return(modelSummary)
  
}

#apply function to an example dataset
fitStanModel(allYrs[which(allYrs$site_id=="100000001"),])

#including year random effects - probably not necessary (recommened by Daskalova et al.)
#maybe check later
#fit1 <- brm(spp_richness ~ year_wMissing + (1|year_wMissing) + ar(time = year_wMissing, p = 1),data = site100000001, family = poisson())


##load pre-calculated slopes
slopeData <- read.csv("outputs/All_siteLevel_and_glmOutput.csv", header=T) # change file name according to the time series to be analyzed
head(slopeData)

#combine trends in mixed model

#Species richness dataset
sr <- slopeData[!is.na(slopeData$SppRich_Est),]
#where w = 1/sd of the trend estimates

sr$SppRich_weights <- 1/sr$SppRich_SE
hist(sr$SppRich_weights)

#look at normality of slopes
hist(sr$SppRich_Est)

#what random effects?
#country?
#site
table(sr$site)
#doesnt make sense to include site as a random ID since we only have one value per site


#see what the default priors are
get_prior(SppRich_Est|weights(SppRich_weights) ~ 1 + (1|country),
          data = sr, family = gaussian())
#default ones ok  - we might play later with this more

prior2 = c(set_prior("cauchy(0,2)", class = "sd"))#cauchy prior common for sd 

fit1 <- brm(SppRich_Est|weights(SppRich_weights) ~ 1 + (1|country),
            data = sr, family = gaussian(), prior = prior2)

summary(fit1)
plot(fit1)

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