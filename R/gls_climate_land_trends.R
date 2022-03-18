##Set working directory
setwd("C:/Users/Ellen/Desktop/aquatic_data")

# load libraries
library(nlme)

# attach data
DATA1_list <- read.csv("All_indices_benthicMacroInverts_AllYears.csv", header=T) # change file name according to the time series to be analyzed
attach(DATA1_list)
head(DATA1_list)
turnover<-as.numeric(turnover)

##############################make more date variables
#centre Year - helps model convergence to center variables for the model
DATA1_list$cYear <- DATA1_list$year - median(DATA1_list$year)

#or just have as an index starting from 1
DATA1_list$iYear <- DATA1_list$year - min(DATA1_list$year)+1

#centre day of year
DATA1_list$cday_of_year <- DATA1_list$day_of_year - median(DATA1_list$day_of_year,na.rm=T)

####################################
options(scipen = 999)
#################################################################
################################################################################################
#CLIMATE############################################
#aet
hist((aet_mm_12moPrior))

df<-data.frame(site_id, iYear, cYear, aet_mm_12moPrior, ppt_mm_12moPrior, tmin_C_12moPrior, tmax_C_12moPrior)
x <- c("site", "iYear", "cYear", "aet", "ppt", "tmin", "tmax")
    colnames(df) <- x
head(df)

trends <- NULL
for(i in unique(df$site)){
  tryCatch({
  sub <- df[df$site == i, ]
  trend.i <- summary(gls(aet ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

#order site results
aet_df <- trends[order(trends$site),]

#rename gls output
xn <- c("aet_site", "aet_Est", "aet_SE", "aet_p")
    colnames(aet_df) <- xn
head(aet_df)
nrow(aet_df)
############################################################
##############################################
#ppt
hist((ppt_mm_12moPrior))

trends <- NULL
for(i in unique(df$site)){
  tryCatch({
  sub <- df[df$site == i, ]
  trend.i <- summary(gls(ppt ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

#order site results
ppt_df <- trends[order(trends$site),]

#rename gls output
xn <- c("ppt_site", "ppt_Est", "ppt_SE", "ppt_p")
    colnames(ppt_df) <- xn
head(ppt_df)
nrow(ppt_df)
##############################################
#tmax
hist((tmax_C_12moPrior))

trends <- NULL
for(i in unique(df$site)){
  tryCatch({
  sub <- df[df$site == i, ]
  trend.i <- summary(gls(tmax ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

#order site results
tmax_df <- trends[order(trends$site),]

#rename gls output
xn <- c("tmax_site", "tmax_Est", "tmax_SE", "tmax_p")
    colnames(tmax_df) <- xn
head(tmax_df)
nrow(tmax_df)
##############################################
#tmin
hist((tmin_C_12moPrior))

trends <- NULL
for(i in unique(df$site)){
  tryCatch({
  sub <- df[df$site == i, ]
  trend.i <- summary(gls(tmin ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

#order site results
tmin_df <- trends[order(trends$site),]

#rename gls output
xn <- c("tmin_site", "tmin_Est", "tmin_SE", "tmin_p")
    colnames(tmin_df) <- xn
head(tmin_df)
nrow(tmin_df)
##############################################
##############################################

climate.tr<- cbind(aet_df, ppt_df,tmax_df,tmin_df)
head(climate.tr)

identical(climate.tr$aet_site,climate.tr$ppt_site)
identical(climate.tr$aet_site,climate.tr$tmax_site)
identical(climate.tr$aet_site,climate.tr$tmin_site)

library(dplyr)
new_df<-select(climate.tr, -c("ppt_site","tmax_site","tmin_site"))

names(new_df)[names(new_df) == 'aet_site'] <- 'site'
head(new_df)

##
write.csv(new_df,"Climate_glmTrends.csv")

#########################################################4
###########################################################################################
################################################################################################
#NITROGEN############################################
#NH4
hist((NH4))

df<-data.frame(site_id, iYear, cYear, NH4, NO3, N)
x <- c("site", "iYear", "cYear", "NH4", "NO3", "N")
    colnames(df) <- x
head(df2)
df2 <- df[!is.na(df$NH4),]

#total values for each site to identify sites with no unique values over time
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = length(unique((sub$NH4)))
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 1),]
nrow(badi)
goodi <- range[ which(range$total.i > 1),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#length of each site to pull out time series with <4 obs
range <- NULL
for(i in unique(df3$site)){
  sub <- df3[df3$site == i, ]
  total.i = nrow(sub)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi2 <- range[ which(range$total.i < 4),]
nrow(badi2)
goodi2 <- range[ which(range$total.i >= 4),]

# Create site subset for sites with at least 4 obs
check = (df3$site %in% goodi2$site) #make a list of T or F is there is a match in spp
df4 = df3[check,] #subset by this list

trends <- NULL
for(i in unique(df4$site)){
  tryCatch({
  sub <- df4[df4$site == i, ]
  trend.i <- summary(gls(NH4 ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

bla <- rep("NA",each=nrow(badi2))
zeries2 <- data.frame(badi2$site, bla, bla, bla) 
colnames(zeries2) <- xn

#order site results
tot_t<-rbind(trends, zeries, zeries2)
NH4_df <- tot_t[order(tot_t$site),]

#rename gls output
xn <- c("NH4_site", "NH4_Est", "NH4_SE", "NH4_p")
    colnames(NH4_df) <- xn
head(NH4_df)
nrow(NH4_df)
############################################################
#NO3
hist((NO3))

#total values for each site to identify sites with no unique values over time
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = length(unique((sub$NO3)))
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 1),]
nrow(badi)
goodi <- range[ which(range$total.i > 1),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#length of each site to pull out time series with <4 obs
range <- NULL
for(i in unique(df3$site)){
  sub <- df3[df3$site == i, ]
  total.i = nrow(sub)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi2 <- range[ which(range$total.i < 4),]
nrow(badi2)
goodi2 <- range[ which(range$total.i >= 4),]

# Create site subset for sites with at least 4 obs
check = (df3$site %in% goodi2$site) #make a list of T or F is there is a match in spp
df4 = df3[check,] #subset by this list

trends <- NULL
for(i in unique(df4$site)){
  tryCatch({
  sub <- df4[df4$site == i, ]
  trend.i <- summary(gls(NO3 ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

bla <- rep("NA",each=nrow(badi2))
zeries2 <- data.frame(badi2$site, bla, bla, bla) 
colnames(zeries2) <- xn

#order site results
tot_t<-rbind(trends, zeries, zeries2)
NO3_df <- tot_t[order(tot_t$site),]

#rename gls output
xn <- c("NO3_site", "NO3_Est", "NO3_SE", "NO3_p")
    colnames(NO3_df) <- xn
head(NO3_df)
nrow(NO3_df)
######################################################################
#N
hist((N))

#total values for each site to identify sites with no unique values over time
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = length(unique((sub$N)))
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 1),]
nrow(badi)
goodi <- range[ which(range$total.i > 1),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#length of each site to pull out time series with <4 obs
range <- NULL
for(i in unique(df3$site)){
  sub <- df3[df3$site == i, ]
  total.i = nrow(sub)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi2 <- range[ which(range$total.i < 4),]
nrow(badi2)
goodi2 <- range[ which(range$total.i >= 4),]

# Create site subset for sites with at least 5 obs
check = (df3$site %in% goodi2$site) #make a list of T or F is there is a match in spp
df4 = df3[check,] #subset by this list

trends <- NULL
for(i in unique(df4$site)){
  tryCatch({
  sub <- df4[df4$site == i, ]
  trend.i <- summary(gls(N ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

bla <- rep("NA",each=nrow(badi2))
zeries2 <- data.frame(badi2$site, bla, bla, bla) 
colnames(zeries2) <- xn

#order site results
tot_t<-rbind(trends, zeries, zeries2)
N_df <- tot_t[order(tot_t$site),]

#rename gls output
xn <- c("N_site", "N_Est", "N_SE", "N_p")
    colnames(N_df) <- xn
head(N_df)
nrow(N_df)
##############################################
##############################################

nitrogen.tr<- cbind(NH4_df, NO3_df,N_df)
head(nitrogen.tr)

identical(nitrogen.tr$NH4_site,nitrogen.tr$NO3_site)
identical(nitrogen.tr$NH4_site,nitrogen.tr$N_site)

library(dplyr)
new_df<-select(nitrogen.tr, -c("NO3_site","N_site"))

names(new_df)[names(new_df) == 'NH4_site'] <- 'site_ID'
head(new_df)

##
write.csv(new_df,"Nitrogen_glmTrends.csv")

##
#######################################################################
#################################################################
################################################################################################
#LAND COVER############################################
#crop
hist((crop_perc_upstr))

df<-data.frame(site_id, iYear, cYear, crop_perc_upstr, forest_perc_upstr, grassShrubland_perc_upstr, urban_perc_upstr, water_perc_upstr)
x <- c("site", "iYear", "cYear", "crop", "forest", "grass", "urban", "water")
    colnames(df) <- x
head(df2)
df2 <- df[!is.na(df$crop),]

#total values for each site to identify sites with no unique values over time
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = length(unique((sub$crop)))
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 1),]
nrow(badi)
goodi <- range[ which(range$total.i > 1),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#length of each site to pull out time series with <4 obs
range <- NULL
for(i in unique(df3$site)){
  sub <- df3[df3$site == i, ]
  total.i = nrow(sub)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi2 <- range[ which(range$total.i < 4),]
nrow(badi2)
goodi2 <- range[ which(range$total.i >= 4),]

# Create site subset for sites with at least 4 obs
check = (df3$site %in% goodi2$site) #make a list of T or F is there is a match in spp
df4 = df3[check,] #subset by this list

trends <- NULL
for(i in unique(df4$site)){
  tryCatch({
  sub <- df4[df4$site == i, ]
  trend.i <- summary(gls(crop ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

bla <- rep("NA",each=nrow(badi2))
zeries2 <- data.frame(badi2$site, bla, bla, bla) 
colnames(zeries2) <- xn

#order site results
tot_t<-rbind(trends, zeries, zeries2)
crop_df <- tot_t[order(tot_t$site),]

#rename gls output
xn <- c("crop_site", "crop_Est", "crop_SE", "crop_p")
    colnames(crop_df) <- xn
head(crop_df)
nrow(crop_df)
##############################################
##################################################################################
##forest
#total values for each site to identify sites with no unique values over time
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = length(unique((sub$forest)))
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 1),]
nrow(badi)
goodi <- range[ which(range$total.i > 1),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#length of each site to pull out time series with <4 obs
range <- NULL
for(i in unique(df3$site)){
  sub <- df3[df3$site == i, ]
  total.i = nrow(sub)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi2 <- range[ which(range$total.i < 4),]
nrow(badi2)
goodi2 <- range[ which(range$total.i >= 4),]

# Create site subset for sites with at least 4 obs
check = (df3$site %in% goodi2$site) #make a list of T or F is there is a match in spp
df4 = df3[check,] #subset by this list

trends <- NULL
for(i in unique(df4$site)){
  tryCatch({
  sub <- df4[df4$site == i, ]
  trend.i <- summary(gls(forest ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

bla <- rep("NA",each=nrow(badi2))
zeries2 <- data.frame(badi2$site, bla, bla, bla) 
colnames(zeries2) <- xn

#order site results
tot_t<-rbind(trends, zeries, zeries2)
forest_df <- tot_t[order(tot_t$site),]

#rename gls output
xn <- c("forest_site", "forest_Est", "forest_SE", "forest_p")
    colnames(forest_df) <- xn
head(forest_df)
nrow(forest_df)
##############################################
##################################################################################
##grassland/shrubland
#total values for each site to identify sites with no unique values over time
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = length(unique((sub$grass)))
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 1),]
nrow(badi)
goodi <- range[ which(range$total.i > 1),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#length of each site to pull out time series with <4 obs
range <- NULL
for(i in unique(df3$site)){
  sub <- df3[df3$site == i, ]
  total.i = nrow(sub)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi2 <- range[ which(range$total.i < 4),]
nrow(badi2)
goodi2 <- range[ which(range$total.i >= 4),]

# Create site subset for sites with at least 4 obs
check = (df3$site %in% goodi2$site) #make a list of T or F is there is a match in spp
df4 = df3[check,] #subset by this list

trends <- NULL
for(i in unique(df4$site)){
  tryCatch({
  sub <- df4[df4$site == i, ]
  trend.i <- summary(gls(grass ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

bla <- rep("NA",each=nrow(badi2))
zeries2 <- data.frame(badi2$site, bla, bla, bla) 
colnames(zeries2) <- xn

#order site results
tot_t<-rbind(trends, zeries, zeries2)
grass_df <- tot_t[order(tot_t$site),]

#rename gls output
xn <- c("grass_site", "grass_Est", "grass_SE", "grass_p")
    colnames(grass_df) <- xn
head(grass_df)
nrow(grass_df)
##############################################
##################################################################################
##urban
#total values for each site to identify sites with no unique values over time
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = length(unique((sub$urban)))
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 1),]
nrow(badi)
goodi <- range[ which(range$total.i > 1),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#length of each site to pull out time series with <4 obs
range <- NULL
for(i in unique(df3$site)){
  sub <- df3[df3$site == i, ]
  total.i = nrow(sub)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi2 <- range[ which(range$total.i < 4),]
nrow(badi2)
goodi2 <- range[ which(range$total.i >= 4),]

# Create site subset for sites with at least 4 obs
check = (df3$site %in% goodi2$site) #make a list of T or F is there is a match in spp
df4 = df3[check,] #subset by this list

trends <- NULL
for(i in unique(df4$site)){
  tryCatch({
  sub <- df4[df4$site == i, ]
  trend.i <- summary(gls(urban ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

bla <- rep("NA",each=nrow(badi2))
zeries2 <- data.frame(badi2$site, bla, bla, bla) 
colnames(zeries2) <- xn

#order site results
tot_t<-rbind(trends, zeries, zeries2)
urban_df <- tot_t[order(tot_t$site),]

#rename gls output
xn <- c("urban_site", "urban_Est", "urban_SE", "urban_p")
    colnames(urban_df) <- xn
head(urban_df)
nrow(urban_df)
##############################################
##################################################################################
##water
#total values for each site to identify sites with no unique values over time
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = length(unique((sub$water)))
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 1),]
nrow(badi)
goodi <- range[ which(range$total.i > 1),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#length of each site to pull out time series with <4 obs
range <- NULL
for(i in unique(df3$site)){
  sub <- df3[df3$site == i, ]
  total.i = nrow(sub)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi2 <- range[ which(range$total.i < 4),]
nrow(badi2)
goodi2 <- range[ which(range$total.i >= 4),]

# Create site subset for sites with at least 4 obs
check = (df3$site %in% goodi2$site) #make a list of T or F is there is a match in spp
df4 = df3[check,] #subset by this list

trends <- NULL
for(i in unique(df4$site)){
  tryCatch({
  sub <- df4[df4$site == i, ]
  trend.i <- summary(gls(water ~ cYear, correlation = corAR1(form = ~ iYear),data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

bla <- rep("NA",each=nrow(badi2))
zeries2 <- data.frame(badi2$site, bla, bla, bla) 
colnames(zeries2) <- xn

#order site results
tot_t<-rbind(trends, zeries, zeries2)
water_df <- tot_t[order(tot_t$site),]

#rename gls output
xn <- c("water_site", "water_Est", "water_SE", "water_p")
    colnames(water_df) <- xn
head(water_df)
nrow(water_df)
##############################################
#########################################################################################################################
##############################################

glm.output<- cbind(crop_df,forest_df,grass_df,urban_df,water_df)
head(glm.output)

##glm.output$EPT_E10_site<-as.numeric(glm.output$EPT_E10_site)

identical(glm.output$crop_site,glm.output$forest_site)
identical(glm.output$crop_site,glm.output$grass_site)
identical(glm.output$crop_site,glm.output$urban_site)
identical(glm.output$crop_site,glm.output$water_site)

library(dplyr)
new_df<-select(glm.output, -c("forest_site","grass_site","urban_site","water_site"))

names(new_df)[names(new_df) == 'crop_site'] <- 'site'
head(new_df)

##
write.csv(new_df,"landcover_glmOutput.csv")
##

################################################################
###########################################################