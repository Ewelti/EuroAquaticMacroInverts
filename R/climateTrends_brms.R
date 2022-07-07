##Set working directory
setwd("C:/Users/Ellen/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

#install.packages("brms")
# if (!require(devtools)) {
#   install.packages("devtools")
# }
#install.packages("rlang")
#install.packages("RcppParallel")
#install.packages("RcppEigen")

library(rstan)
library(brms)
library(lubridate)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#load data
allYrs <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T)
head(allYrs)

##############################make more date variables
#centre Year - helps model convergence to center variables for the model
allYrs$cYear <- allYrs$year - median(allYrs$year)

#or just have as an index starting from 1
allYrs$iYear <- allYrs$year - min(allYrs$year)+1
###############################################################

hist(allYrs$ppt_mm_12moPrior)
hist(allYrs$tmax_C_12moPrior)

#### two-stage models ####

## check one site
#fit trend model to a single time-series
site100000001<-allYrs[which(allYrs$site_id=="100000001"),]

myD <- make_standata(ppt_mm_12moPrior ~ cYear, data = site100000001)
myD$meanResponse <- round(median(site100000001$ppt_mm_12moPrior), 1)
myD$sdResponse <- max(round(mad(site100000001$ppt_mm_12moPrior), 1), 2.5)
fit1 <-stan('climate_stan_model.stan', 
            data = myD, chains = 4,iter = 1000)
trend.i <- (extract(fit1, pars="b[1]",permuted = FALSE))
trend.i <- as.matrix(fit1, pars = c("b[1]"))
head(trend.i)
fit_summary <- summary(fit1, pars = c("b[1]"))$summary
trend.i <- fit_summary[,1]

##############################################

##### Precipitation trends

trends <- NULL

for(i in unique(allYrs$site_id)){
  
  #subset dataset to focal site
  sub <- allYrs[allYrs$site_id == i, ]
  
  #remove NAs
  sub <- subset(sub, !is.na(ppt_mm_12moPrior))
  
  #define model -only needed on first run - after that we call the file
  #prior1 = c(set_prior("normal(0,10)", class = "b"))
  #myCode <- make_stancode(ppt_mm_12moPrior ~ cYear, data = sub, 
  #                        prior = prior1, save_model = "climate_stan_model")
  
  #define data
  myData <- make_standata(ppt_mm_12moPrior ~ cYear, data = sub)
  myData$meanResponse <- round(median(sub$ppt_mm_12moPrior), 1)
  myData$sdResponse <- max(round(mad(sub$ppt_mm_12moPrior), 1), 2.5)
  
  #run model
  fit1 <-stan('outputs/climate_stan_model.stan', 
              data = myData, chains = 4,iter = 1000)
  
  fit_summary <- summary(fit1, pars = c("b[1]"))$summary
  trend.i <- fit_summary[,1]
  trend.i <- data.frame(site_id = i, trend = trend.i)
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
  
}

#order site results
ppt_df <- trends[order(trends$site_id),]
head(ppt_df)

#rename output variables
xn <- c("site", "ppt_Est")
colnames(ppt_df) <- xn
nrow(ppt_df)

write.csv(ppt_df, "precipitation_trends.csv")

### temperature

trends <- NULL

for(i in unique(allYrs$site_id)){
  
  #subset dataset to focal site
  sub <- allYrs[allYrs$site_id == i, ]
  
  #remove NAs
  sub <- subset(sub, !is.na(tmax_C_12moPrior))
  
  #define model -only needed on first run - after that we call the file
  #prior1 = c(set_prior("normal(0,1)", class = "b"))
  #myCode <- make_stancode(tmax_C_12moPrior ~ cYear, data = sub, 
  #                        prior = prior1, save_model = "climate_stan_model")
  
  #define data
  myData <- make_standata(tmax_C_12moPrior ~ cYear, data = sub)
  myData$meanResponse <- round(median(sub$tmax_C_12moPrior), 1)
  myData$sdResponse <- max(round(mad(sub$tmax_C_12moPrior), 1), 2.5)
  
  #run model
  fit1 <-stan('outputs/climate_stan_model.stan', 
              data = myData, chains = 4,iter = 1000)
  
  fit_summary <- summary(fit1, pars = c("b[1]"))$summary
  trend.i <- fit_summary[,1]
  trend.i <- data.frame(site_id = i, trend = trend.i)
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
  
}

#order site results
tmax_df <- trends[order(trends$site_id),]
head(tmax_df)

#rename output variables
xn <- c("site", "tmax_Est")
colnames(tmax_df) <- xn
nrow(tmax_df)

write.csv(tmax_df, "temperature_trends.csv")