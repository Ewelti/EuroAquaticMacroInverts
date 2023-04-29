##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

# attach data
response_stan <- readRDS("outputs/outputs_metaAnalysis/stanTrends_site_level.rds")
head(response_stan)
unique(response_stan$Response)
length(unique(response_stan$site_id))

metricdata <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv")
head(metricdata)

one.ESTs <- read.csv("outputs/outputs_sensitivity/metaanalysisModelComparison/Yr_OneStage_Ests.csv")
head(one.ESTs)

F_to_Est <- subset(one.ESTs, Response == "func_turnover")
head(F_to_Est)

##pull out min and max years and F_to in yr 1 as intercept estimate
FT_Data <- unique(metricdata[,c("site_id","F_to", "year_wMissing")])
FT_Data <- na.omit(FT_Data)

library(data.table)
FT_Datam <- data.table(FT_Data)
####mins <- FT_Datam[ , .SD[which.min(year_wMissing)], by = site_id]
mins <- FT_Datam[ , list(year_wMissing = min(year_wMissing)), by = site_id]
colnames(mins)[colnames(mins) == "year_wMissing"] <- "min_yr"
max <- FT_Datam[ , list(year_wMissing = max(year_wMissing)), by = site_id]
colnames(max)[colnames(max) == "year_wMissing"] <- "max_yr"
ft_me <- merge(mins,max,by="site_id")
med <- FT_Datam[ , list(year_wMissing = median(year_wMissing)), by = site_id]
colnames(med)[colnames(med) == "year_wMissing"] <- "med_yr"
ft_mer <- merge(ft_me,med,by="site_id")
inter <- FT_Datam[ , list(F_to = mean(F_to)), by = site_id]
colnames(inter)[colnames(inter) == "F_to"] <- "intercept"
ft_merg <- merge(ft_mer,inter,by="site_id")
head(ft_merg)

##predict % change/yr for each sites
library(betareg)

allinter <- mean(ft_merg$intercept)
minyr <- mean(ft_merg$min_yr)
maxyr <- mean(ft_merg$max_yr)

Year.num <- maxyr - minyr
yr.1 <- 0-Year.num/2
yr.n <- Year.num/2

(year1 <- plogis(allinter + F_to_Est$Estimate * yr.1))
(yearn <- plogis(allinter + F_to_Est$Estimate * yr.n))

(rate <- (yearn/year1)^(1/(maxyr-minyr)))
(perc = (rate -1)*100) 

percC <- NULL
for(i in F_to_Est[2:11]){
	year1 <- plogis(allinter + i * yr.1)
	yearn <- plogis(allinter + i * yr.n)
	rate <- (yearn/year1)^(1/(maxyr-minyr))
	perc = (rate -1)*100
    percC.i <- data.frame(perc)
    percC <- rbind(percC, percC.i) ; rm(percC.i, year1, yearn, rate, perc)
} ; rm(i)
percC
stan_Fto <- cbind(colnames(F_to_Est[2:11]),percC)
colnames(stan_Fto)[colnames(stan_Fto) == "colnames(F_to_Est[2:11])"] <- "output"

##############################################################################################
##unweighted
uw.ESTs <- read.csv("outputs/outputs_sensitivity/metaanalysisModelComparison/Yr_metaanaly_Ests_unweighted.csv")
head(uw.ESTs)

F_to_Est <- subset(uw.ESTs, Response == "func_turnover")
head(F_to_Est)

(year1 <- plogis(allinter + F_to_Est$Estimate * yr.1))
(yearn <- plogis(allinter + F_to_Est$Estimate * yr.n))

(rate <- (yearn/year1)^(1/(maxyr-minyr)))
(perc = (rate -1)*100) 

percC <- NULL
for(i in F_to_Est[2:11]){
	year1 <- plogis(allinter + i * yr.1)
	yearn <- plogis(allinter + i * yr.n)
	rate <- (yearn/year1)^(1/(maxyr-minyr))
	perc = (rate -1)*100
    percC.i <- data.frame(perc)
    percC <- rbind(percC, percC.i) ; rm(percC.i, year1, yearn, rate, perc)
} ; rm(i)
percC
stan_Fto <- cbind(colnames(F_to_Est[2:11]),percC)
colnames(stan_Fto)[colnames(stan_Fto) == "colnames(F_to_Est[2:11])"] <- "output"

