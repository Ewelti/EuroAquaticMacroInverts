#install.packages("brms")
# if (!require(devtools)) {
#   install.packages("devtools")
# }


library(brms)
library(lubridate)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#load data
d1 <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears.csv", header=T) # change file name according to the time series to be analyzed
allYrs <- d1[!is.na(d1$site_id_wMissing),]
allYrs$turnover<-as.numeric(allYrs$turnover) # make numeric
head(allYrs)

#centre Year - helps model convergence to center variables for the model
allYrs$cYear <- allYrs$year_wMissing - median(allYrs$year_wMissing)
summary(allYrs$cYear)
#or just have as an index starting from 1
allYrs$iYear <- allYrs$year_wMissing - min(allYrs$year_wMissing)+1
summary(allYrs$iYear)

#centre day of year
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
prior1 = c(set_prior("normal(0,0.5)", class = "ar"),
           set_prior("normal(0,10)", class = "b"))

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

# write a function to consider day of year if sampling is more than 30 days apart
fitStanModel <- function(mydata){
  
  #if sampling occurs in more than one month include a seasonal term in the model
  maxDiffDays = max(mydata$cday_of_year)-min(mydata$cday_of_year)
  
  if(maxDiffDays < 30) {
    myformula <- bf(spp_richness ~ cYear + ar(time = iYear, p = 1))
  } else{
    myformula <- bf(spp_richness ~ cday_of_year + cYear + ar(time = iYear, p = 1))
  }
  
  #fit model
  fit1 <- brm(myformula, data = mydata, family = poisson(), prior = prior1, refresh = 0, cores=1)
  
  #extract model fits
  modelSummary <- fixef(fit1, pars="cYear")[1, c(1,2)]
  return(modelSummary)
  
}

#apply function to an example dataset
est <- fitStanModel(allYrs[which(allYrs$site_id=="100000001"),])

#including year random effects - probably not necessary (recommened by Daskalova et al.)
#maybe check later
#fit1 <- brm(spp_richness ~ year_wMissing + (1|year_wMissing) + ar(time = year_wMissing, p = 1),data = site100000001, family = poisson())

#loop for all sites ####has problems with crashing-- trying tryCatch
trends <- NULL
for(i in unique(allYrs$site_id)){
  tryCatch({
    sub <- allYrs[allYrs$site_id == i, ]
    trend.i <- fitStanModel(sub)
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})    
} ; rm(i)


#write to a text file during model fitting - to see where it is going wrong

for(i in unique(allYrs$site_id)){
  sub <- allYrs[allYrs$site_id == i, ]
  trend.i <- fitStanModel(sub)
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  write.table(trend.i, file="outputs/running_modelOutput.txt", 
              append=TRUE, row.names=FALSE, sep=",") 
}

#check where it crashes
fitStanModel(allYrs[which(allYrs$site_id == unique(allYrs$site_id)[5]),])

##load pre-calculated slopes
slopeData <- read.csv("outputs/All_siteLevel_and_glmOutput.csv", header=T)
head(slopeData)

######################################################################################
##variable for discussion
################################
####possible random effects:
# study_id, Country, 
##redundant with study_id but maybe put as fixed at a coarse level?: sampling_method_long
##too few levels: season, TaxonomicRes
##manybe too many levels: MacrobasinID
################################
####possible fixed effects:

##climate slopes and means:
# aet_Est, ppt_Est, tmax_Est, tmin_Est, aet_mm_12moPrior, ppt_mm_12moPrior, tmax_C_12moPrior, tmin_C_12moPrior

##N slopes and means:
# NH4_Est, NO3_Est, NH4_mean, NO3_mean

##stream characteristics:
# strahler_streamOrder, accumulation_atPoint, elevation_atPoint, slope_atPoint, 

##land use (right now just microbasin):
# crop_meanPerc, forest_meanPerc, grassShrubland_meanPerc, urban_meanPerc

##dams
# dam_impact_score_lessthan100km, dam_minDist_km_lessthan100km, dam_num_connected_lessthan100km # probably only need one of these


#for all fixed effects - first scale the continuous variables by their standard deviations

#function to add a new column onto the data with scaled vars (with s before their name)
scaleVars <- function(df){
  newd <- plyr::numcolwise(scale)(df)
  names(newd) <- sapply(names(newd),function(x)paste0("s",x))
  cbind(df, newd)
}
#apply function
slopeData <- scaleVars(slopeData)

#Species richness dataset
sr <- slopeData[!is.na(slopeData$SppRich_Est),]
#where w = 1/sd of the trend estimates

sr$SppRich_weights <- 1/sr$SppRich_SE
hist(sr$SppRich_weights)

#look at normality of slopes
hist(sr$SppRich_Est)

#example driver model
prior2 = c(set_prior("cauchy(0,2)", class = "sd"),
           set_prior("normal(0,5)", class = "b"))

fit1 <- brm(SppRich_Est|weights(SppRich_weights) ~ surban_meanPerc + (1|Country) + (1|study_id),
            data = sr, family = gaussian(), prior = prior2)

summary(fit1)

##############################
####possible response variables (these also all have calcuated standard error, coded as 'Response_SE'):
##Taxonomic diversity:
# SppRich_Est, Adun_Est, TurnO_Est
##maybe will not use or swop for other indices:
# SimpD_Est, ShanH_Est, EvenJ_Est

##Functional diversity:
# F_to_Est, FRic_Est, FEve_Est, FDiv_Est, RaoQ_Est

##ALiens (possible supplement, this can only be done on a subset of the data IDed to finer taxonomic level)
# AlienAbun_Est, AlienSppRich_Est, nativeAbun_Est, nativeSppRich_Est
#######################################################################################

#see what the default priors are
get_prior(SppRich_Est|weights(SppRich_weights) ~ 1 + (1|Country) + (1|study_id),
          data = sr, family = gaussian())
#default ones ok  - we might play later with this more

prior2 = c(set_prior("cauchy(0,2)", class = "sd"))#cauchy prior common for sd 

fit1 <- brm(SppRich_Est|weights(SppRich_weights) ~ 1 + (1|study_id), 
            data = sr, family = gaussian(), prior = prior2)
#+ (1|Country) # i think we don't need country as no study spans multiple countries (or could nest study_id in country)

summary(fit1)
plot(fit1)


### including spatial autocorrelation
#?
#We should reproject the long/lat into utm coordinates eventually
# the datum for all should be WGS 84 decimal

#get distance matrix for CAR (spatial conditional autoregressive)
distance <- as.matrix(dist(sr[,c("Longitude_X","Latitude_Y")]))
K <- nrow(sr)
W <- array(0, c(K, K))
W[distance == 1] <- 1 	

#lets fiddle with priors later
fit1 <- brm(SppRich_Est|weights(SppRich_weights) ~ 1 + (1|study_id) + car(W, type="icar"),
            data = sr, data2 = list(W = W), 
            family = gaussian())


### one-stage model #####

prior1 = c(set_prior("normal(0,0.5)", class = "ar"),
           set_prior("normal(0,5)", class = "b"),
           set_prior("cauchy(0,2)", class = "sd"))


fit1 <- brm(spp_richness ~ cYear + (1 + cYear|study_id) + 
              ar(time = iYear, p = 1), 
            data = allYrs, family = poisson(),
            niter = 100)


### moving average ####

fitMAtrend <- function(allYrs, startYear, timespan = 10){
  
  #restrict to time period of interest
  allYrsS <- subset(allYrs, year_wMissing >= startYear)
  allYrsS <- subset(allYrsS, year_wMissing < (startYear+timespan))
  
  #restrict to study with sufficient data in this time period - 5 years?
  study_periods <- tapply(allYrsS$year_wMissing,allYrsS$site_id,
                          function(x)length(unique(x)))
  allYrsS <- subset(allYrsS, site_id %in% 
                      names(study_periods)[study_periods>=5])                 
  #fit model to this subset
  fit1 <- brm(spp_richness ~ cYear, 
              #+ (1 + cYear|study_id) 
              #+ ar(time = iYear, p = 1),
              data = allYrsS, 
              family = poisson())
  
  #extract trend coefficient
  temp <- data.frame(fixef(fit1, pars="cYear"))
  temp$startYear = startYear
  temp$timespan = timespan
  return(temp)
  
}

#run model

#decide time span to run model over
table(allYrs$year_wMissing)
years <- 1980:2010

#run model over these years
modelsList <- lapply(years, function(i){
  fitMAtrend(allYrs, startYear = i, timespan = 10)
})
maTrends <- do.call(rbind,modelsList)


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