##Set working directory
#setwd("C:/Users/ewelti/Desktop/aquatic_data")

#install.packages("brms")
# if (!require(devtools)) {
#   install.packages("devtools")
# }
#install.packages("rlang")

#### Attempt to find out what is crashing brms ####
## using spp richness and those sites with sampling within a 30 day period as an example
## unfortunately crashes seem random- crashing on different sites each run :/


library(brms)
library(lubridate)

#load data
d1 <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears.csv", header=T) # change file name according to the time series to be analyzed
allYrs <- d1[!is.na(d1$site_id_wMissing),]
allYrs$estsnover<-as.numeric(allYrs$estsnover) # make numeric
head(allYrs)

#centre Year - helps model convergence to center variables for the model
allYrs$cYear <- allYrs$year_wMissing - median(allYrs$year_wMissing)

#or just have as an index starting from 1
allYrs$iYear <- allYrs$year_wMissing - min(allYrs$year_wMissing)+1

#centre day of year
allYrs$cday_of_year <- allYrs$day_of_year - median(allYrs$day_of_year,na.rm=T)
  
#### two-stage models ####

#set priors now
prior1 = c(set_prior("normal(0,0.5)", class = "ar"), 
          set_prior("normal(0,5)", class = "b"))

############################################################
##break up sites by those with sampling is more than 30 days apart
#calculate range of day of year for each site

#pull out range of days of year
range <- NULL
for(i in unique(allYrs$site_id)){
  sub <- allYrs[allYrs$site_id == i, ]
  maxDiffDays = max(sub$cday_of_year)-min(sub$cday_of_year)
  maxDiffDays <- data.frame(site = i, 
                        t(maxDiffDays))
  range <- rbind(range, maxDiffDays) ; rm(maxDiffDays, sub)
} ; rm(i)

#which sites have sampling less and more than 30 days apart
LT30 <- range[ which(range$t.maxDiffDays. <= 30),]
MT30 <- range[ which(range$t.maxDiffDays. > 30),]

# subset all data for sites with with sampling ranging less than 30 days apart
allyrs_LT30 = (allYrs$site_id %in% LT30$site) #make a list of T or F is there is a match in spp
ay_LT30 = allYrs[allyrs_LT30,] #subset by this list
head(ay_LT30)
unique(ay_LT30$site_id)
LT30_sppR_un <- unique(ay_LT30$site_id)
length(LT30_sppR_un)

hist(ay_LT30$spp_richness)
LT30_sppR <- ay_LT30[!is.na(ay_LT30$spp_richness),]

myformula <- bf(spp_richness ~ cYear + ar(time = iYear, p = 1))


site <- LT30_sppR_un[1]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests_b<-cbind(site[1],trend.i[,1],trend.i[,2])
colnames(ests_b) <- c("site_id", "SppRich_Est", "SppRich_Est.Err")

site <- LT30_sppR_un[2]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b)
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[3]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[4]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[5]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[6]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

#site <- LT30_sppR_un[7]
#sub <- LT30_sppR[LT30_sppR$site_id == site, ]
#fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
#trend.i <- fixef(fit1, pars="cYear")
#ests<-cbind(site[1],trend.i[,1],trend.i[,2])
#ests_b<-rbind(ests,ests_b) 
#write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[8]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[9]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[10]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[11]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[12]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[13]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[14]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[15]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[16]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[17]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[18]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[19]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[20]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[21]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[22]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[23]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[24]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[25]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[26]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[27]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[28]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[29]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[30]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[31]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[32]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[33]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[34]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[35]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[36]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[37]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[38]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[39]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[40]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[41]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[42]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[43]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[44]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[45]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[46]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[47]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[48]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[49]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[50]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[51]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[52]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[53]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[54]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[55]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[56]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[57]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[58]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[59]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[60]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[61]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[62]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[63]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[64]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[65]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[66]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[67]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[68]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[69]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[70]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[71]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[72]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[73]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[74]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[75]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[76]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[77]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[78]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[79]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[80]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[81]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[82]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[83]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[84]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[85]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[86]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[87]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[88]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[89]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[90]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[91]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[92]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[93]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[94]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[95]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[96]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[97]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[98]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[99]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[100]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[101]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests_b<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[102]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[103]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[104]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[105]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[106]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[107]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[108]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[109]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[110]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[111]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[112]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[113]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[114]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[115]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[116]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[117]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[118]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[119]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[120]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[121]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[122]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[123]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[124]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[125]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[126]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[127]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[128]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[129]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[130]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[131]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[132]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

site <- LT30_sppR_un[133]
sub <- LT30_sppR[LT30_sppR$site_id == site, ]
fit1 <- brm(myformula, data = sub, family = poisson(), prior = prior1, refresh = 0)
trend.i <- fixef(fit1, pars="cYear")
ests<-cbind(site[1],trend.i[,1],trend.i[,2])
ests_b<-rbind(ests,ests_b) 
write.csv(ests_b,"LT30_sppR_Est.csv")

######################################################################################
##
