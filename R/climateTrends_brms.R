##Set working directory
setwd("C:/Users/Ellen/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs")

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
allYrs <- read.csv("All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T)
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

#simplest model:
fit1 <- brm(ppt_mm_12moPrior ~ cYear, data = site100000001,family = gaussian())
modelSummary <- fixef(fit1, pars="cYear")[1, c(1,2)]

#see what the default priors are
get_prior(ppt_mm_12moPrior ~ cYear, data = site100000001, family = gaussian())

#set priors now
prior1 = c(set_prior("normal(0,0.5)", class = "ar"), #I don't know what/if to change here
          set_prior("normal(0,5)", class = "b"))

#including autocorrelation (of the residuals)
fit1 <- brm(ppt_mm_12moPrior ~ cYear + ar(time = iYear, p = 1),data = site100000001, family = gaussian())
summary(fit1)
modelSummary <- fixef(fit1, pars="cYear")[1, c(1,2)]


#######################################
#### loop to run precip trends for all sites
trends <- NULL
for(i in unique(allYrs$site_id)){
  tryCatch({
  sub <- allYrs[allYrs$site_id == i, ]
  fit1 <- brm(ppt_mm_12moPrior ~ cYear, data = sub, family = gaussian())
  trend.i <- fixef(fit1, pars="cYear")[1, c(1,2)]
  trend.i <- data.frame(site_id = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site_id),conditionMessage(e), "\n")})
} ; rm(i)

#order site results
tmax_df <- trends[order(trends$site_id),]

#rename gls output
xn <- c("ppt_site", "ppt_Est", "ppt_Error")
    colnames(tmax_df) <- xn
head(tmax_df)
nrow(tmax_df)
##############################################

#dianas attempt

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
    fit1 <-stan('climate_stan_model.stan', 
                       data = myData, chains = 4,iter = 1000)
    
    trend.i <- fixef(fit1, pars="cYear")[1,1]
    trend.i <- data.frame(site_id = i, trend = trend.i)
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    
}



######################################################

# try using a function to run model
fitStanModel <- function(mydata){
    myformula <- bf(ppt_mm_12moPrior ~ cYear + ar(time = iYear, p = 1))
  
  	#fit model
  	fit1 <- brm(myformula, data = mydata, family = gaussian(), prior = prior1, refresh = 0)
  
  		#extract model fits
  		modelSummary <- fixef(fit1, pars="cYear")[1, c(1,2)]
  		return(modelSummary)
}

#apply function to an example dataset
e <- fitStanModel(allYrs[which(allYrs$site_id=="100000001"),])


#loop that always crashes
for(i in unique(allYrs$site_id)){
  tryCatch({
    sub <- allYrs[allYrs$site_id == i, ]
    trend.i <- fitStanModel(sub)
    trend.i <- data.frame(site = i, 
                        t(trend.i))
   write.table(trend.i, file="precip_brms_trends.txt", 
              append=TRUE, row.names=FALSE, sep=",") 
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

#write to a text file during model fitting - to see where it is going wrong

for(i in unique(allYrs$site_id)){
  sub <- allYrs[allYrs$site_id == i, ]
  trend.i <- fitStanModel(sub)
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  write.table(trend.i, file="sppRich_modelOutput.txt", 
              append=TRUE, row.names=FALSE, sep=",") 
}
############################################################

