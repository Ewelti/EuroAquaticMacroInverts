##Set working directory
setwd("C:/Users/ewelti/Desktop/aquatic_data")

# load libraries
library(nlme)

# attach data
DATA1_list <- read.csv("All_indices_benthicMacroInverts_AllYears.csv", header=T)
DATA2 <- DATA1_list[!is.na(DATA1_list$site_id_wMissing),]
DATA2$turnover<-as.numeric(DATA2$turnover) # make numeric
head(DATA2)

##############################make more date variables
#centre Year - helps model convergence to center variables for the model
DATA2$cYear <- DATA2$year_wMissing - median(DATA2$year_wMissing)

#or just have as an index starting from 1
DATA2$iYear <- DATA2$year_wMissing - min(DATA2$year_wMissing)+1

#centre day of year
DATA2$cday_of_year <- DATA2$day_of_year - median(DATA2$day_of_year,na.rm=T)

####################################
# 6 datasets have no EPT in whole time series and should be removed for ept: 
#103000121, 103000604, 109000044, 109000049, 109000254, 114000075
ept <- DATA2[which(DATA2$site_id!=103000121& DATA2$site_id!=103000604& DATA2$site_id!=109000044&
DATA2$site_id!=109000049& DATA2$site_id!=109000254& DATA2$site_id!=114000075),]

attach(ept)
head(ept)
#####################################################################
#no sci notation
options(scipen = 999)
#################################################################
#ept spp_richness trend
hist(ept$EPT_SppRich)

#subset data for loop
df<-data.frame(ept$site_id_wMissing, ept$cday_of_year, ept$iYear, ept$cYear, ept$EPT_SppRich)
x <- c("site", "cDOY", "iYear", "cYear", "sr")
    colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(sr ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(sr ~ cDOY + cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[3, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)
head(trends)

df121000087 <- df[which(df$site==121000087),]
df121000087
df121000087 <- data.frame(site = 121000087, Value = 0, Std.Error = "NA", p.value = "NA") #NAs for trends from this site

#order site results
tot_t<-rbind(trends, df121000087)
sr_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("EPT_SppRich_site", "EPT_SppRich_Est", "EPT_SppRich_SE", "EPT_SppRich_p")
    colnames(sr_df) <- xn
head(sr_df)
nrow(sr_df)
mean(sr_df$EPT_SppRich_Est)
################################################
################################################################################################
#ept abundance trend
hist(log10(ept$EPT_Abund))

#subset data for loop
df<-data.frame(ept$site_id_wMissing, ept$cday_of_year, ept$iYear, ept$cYear, log10(ept$EPT_Abund+1))
x <- c("site", "cDOY", "iYear", "cYear", "ab")
    colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(ab ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(ab ~ cDOY + cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[3, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)
head(trends)

df100000295 <- df[which(df$site==100000295),]
max(df100000295$cDOY)-min(df100000295$cDOY)
g1<-gls(ab ~ cYear + cDOY, na.action=na.omit, data = df100000295)
trend.100000295 <- summary(gls(ab ~ cYear + cDOY, na.action=na.omit, data = df100000295))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df100000295 <- data.frame(site = 100000295, t(trend.100000295))

#order site results
tot_t<-rbind(trends, df100000295)
ab_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("EPT_Abund_site", "EPT_Abund_Est", "EPT_Abund_SE", "EPT_Abund_p")
    colnames(ab_df) <- xn
head(ab_df)
nrow(ab_df)

###################################################################################################################################
#ept shannon trend
hist(ept$EPT_Shan)

#subset data for loop
df<-data.frame(ept$site_id_wMissing, ept$cday_of_year, ept$iYear, ept$cYear, ept$EPT_Shan)
x <- c("site", "cDOY", "iYear", "cYear", "sh")
    colnames(df) <- x
head(df)

#total values for each site to identify zeros
range <- NULL
for(i in unique(df$site)){
  sub <- df[df$site == i, ]
  total.i = sum(sub$sh)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 0),]
goodi <- range[ which(range$total.i > 0),]

# Create site subset for sites without all zeros
check = (df$site %in% goodi$site) #make a list of T or F is there is a match in spp
df2 = df[check,] #subset by this list

#calculate gls trends
trends <- NULL
for(i in unique(df2$site)){
  tryCatch({
    sub <- df2[df2$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(sh ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(sh ~ cDOY + cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[3, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)
head(trends)

Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

df100000120 <- df[which(df$site==100000120),]
max(df100000120$cDOY)-min(df100000120$cDOY)
g1<-gls(sh ~ cYear + cDOY, na.action=na.omit, data = df100000120)
trend.100000120 <- summary(gls(sh ~ cYear + cDOY, na.action=na.omit, data = df100000120))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df100000120 <- data.frame(site = 100000120, t(trend.100000120))

df100000241 <- df[which(df$site==100000241),]
max(df100000241$cDOY)-min(df100000241$cDOY)
g1<-gls(sh ~ cYear + cDOY, na.action=na.omit, data = df100000241)
trend.100000241 <- summary(gls(sh ~ cYear + cDOY, na.action=na.omit, data = df100000241))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df100000241 <- data.frame(site = 100000241, t(trend.100000241))

df111000001 <- df[which(df$site==111000001),]
max(df111000001$cDOY)-min(df111000001$cDOY)
g1<-gls(sh ~ cYear, na.action=na.omit, data = df111000001)
trend.111000001 <- summary(gls(sh ~ cYear, na.action=na.omit, data = df111000001))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df111000001 <- data.frame(site = 111000001, t(trend.111000001))

#order site results
tot_t<-rbind(trends, zeries, df100000120, df100000241, df111000001)
sh_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("EPT_Shan_site", "EPT_Shan_Est", "EPT_Shan_SE", "EPT_Shan_p")
    colnames(sh_df) <- xn
head(sh_df)
nrow(sh_df)
#########################################################################################################################
###################################################################################################################################
#ept E10 trend
hist(ept$EPT_E10)
#subset data for loop
df<-data.frame(ept$site_id_wMissing, ept$cday_of_year, ept$iYear, ept$cYear, ept$EPT_E10)
x <- c("site", "cDOY", "iYear", "cYear", "e10")
    colnames(df) <- x
head(df)
df2 <- df[!is.na(df$e10),]

#total values for each site to identify zeros
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = sum(sub$e10)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 0),]
goodi <- range[ which(range$total.i > 0),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#length of each site to pull out time series with <5 obs
range <- NULL
for(i in unique(df3$site)){
  sub <- df3[df3$site == i, ]
  total.i = nrow(sub)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi2 <- range[ which(range$total.i < 5),]
goodi2 <- range[ which(range$total.i >= 5),]

# Create site subset for sites with at least 5 obs
check = (df3$site %in% goodi2$site) #make a list of T or F is there is a match in spp
df4 = df3[check,] #subset by this list

#calculate gls trends
trends <- NULL
for(i in unique(df4$site)){
  tryCatch({
    sub <- df4[df4$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(e10 ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(e10 ~ cDOY + cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[3, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)
head(trends)

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
e10_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("EPT_E10_site", "EPT_E10_Est", "EPT_E10_SE", "EPT_E10_p")
    colnames(e10_df) <- xn
head(e10_df)
nrow(e10_df)
mean(e10_df$EPT_E10_Est)
#########################################################################################################################
###################################################################################################################################
#ept Turnover trend
hist(ept$EPT_Turnover)
#subset data for loop
df<-data.frame(ept$site_id_wMissing, ept$cday_of_year, ept$iYear, ept$cYear, ept$EPT_Turnover)
x <- c("site", "cDOY", "iYear", "cYear", "Turnover")
    colnames(df) <- x
head(df)
df2 <- df[!is.na(df$Turnover),]

#check if any sites where removed
remsites <- setdiff(unique(df$site),unique(df2$site))

#total values for each site to identify zeros
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = length(unique((sub$Turnover)))
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 1),]
goodi <- range[ which(range$total.i > 1),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#length of each site to pull out time series with <5 obs
range <- NULL
for(i in unique(df3$site)){
  sub <- df3[df3$site == i, ]
  total.i = nrow(sub)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi2 <- range[ which(range$total.i < 5),]
goodi2 <- range[ which(range$total.i >= 5),]

# Create site subset for sites with at least 5 obs
check = (df3$site %in% goodi2$site) #make a list of T or F is there is a match in spp
df4 = df3[check,] #subset by this list

#calculate gls trends
trends <- NULL
for(i in unique(df4$site)){
  tryCatch({
    sub <- df4[df4$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(Turnover ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(Turnover ~ cDOY + cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[3, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)
head(trends)
##

Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

bla <- rep("NA",each=nrow(badi2))
zeries2 <- data.frame(badi2$site, bla, bla, bla) 
colnames(zeries2) <- xn

bla <- rep("NA",each=length(remsites))
remies <- data.frame(remsites, bla, bla, bla) 
colnames(remies) <- xn

df107000237 <- df[which(df$site==107000237),]
max(df107000237$cDOY)-min(df107000237$cDOY)
g1<-gls(Turnover ~ cYear + cDOY, na.action=na.omit, data = df107000237)
trend.107000237 <- summary(gls(Turnover ~ cYear+ cDOY, na.action=na.omit, data = df107000237))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df107000237 <- data.frame(site = 107000237, t(trend.107000237))

df108000021 <- df[which(df$site==108000021),]
max(df108000021$cDOY)-min(df108000021$cDOY)
g1<-gls(Turnover ~ cYear + cDOY, na.action=na.omit, data = df108000021)
trend.108000021 <- summary(gls(Turnover ~ cYear + cDOY, na.action=na.omit, data = df108000021))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000021 <- data.frame(site = 108000021, t(trend.108000021))

df108000130 <- df[which(df$site==108000130),]
max(df108000130$cDOY)-min(df108000130$cDOY)
g1<-gls(Turnover ~ cYear + cDOY, na.action=na.omit, data = df108000130)
trend.108000130 <- summary(gls(Turnover ~ cYear + cDOY, na.action=na.omit, data = df108000130))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000130 <- data.frame(site = 108000130, t(trend.108000130))

#order site results
tot_t<-rbind(trends, remies, zeries, zeries2, df107000237, df108000021, df108000130)
to_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("EPT_Turnover_site", "EPT_Turnover_Est", "EPT_Turnover_SE", "EPT_Turnover_p")
    colnames(to_df) <- xn
head(to_df)
nrow(to_df)

###################################################################################################################################
#ept S_PIE trend
hist(log10(ept$EPT_S_PIE))
#subset data for loop
df<-data.frame(ept$site_id_wMissing, ept$cday_of_year, ept$iYear, ept$cYear, log10(ept$EPT_S_PIE+1))
x <- c("site", "cDOY", "iYear", "cYear", "S_PIE")
    colnames(df) <- x
head(df)
df2 <- df[!is.na(df$S_PIE),]

#check if any sites where removed
remsites <- setdiff(unique(df$site),unique(df2$site))

#total values for each site to identify zeros
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = length(unique((sub$S_PIE)))
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 1),]
goodi <- range[ which(range$total.i > 1),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#length of each site to pull out time series with <5 obs
range <- NULL
for(i in unique(df3$site)){
  sub <- df3[df3$site == i, ]
  total.i = nrow(sub)
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi2 <- range[ which(range$total.i < 5),]
goodi2 <- range[ which(range$total.i >= 5),]

# Create site subset for sites with at least 5 obs
check = (df3$site %in% goodi2$site) #make a list of T or F is there is a match in spp
df4 = df3[check,] #subset by this list

#calculate gls trends
trends <- NULL
for(i in unique(df4$site)){
  tryCatch({
    sub <- df4[df4$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(S_PIE ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(S_PIE ~ cDOY + cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[3, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)
head(trends)
##
 
Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

bla <- rep("NA",each=nrow(badi2))
zeries2 <- data.frame(badi2$site, bla, bla, bla) 
colnames(zeries2) <- xn

bla <- rep("NA",each=length(remsites))
remies <- data.frame(remsites, bla, bla, bla) 
colnames(remies) <- xn

df100000034 <- df[which(df$site==100000034),]
max(df100000034$cDOY)-min(df100000034$cDOY)
g1<-gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df100000034)
trend.100000034 <- summary(gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df100000034))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df100000034 <- data.frame(site = 100000034, t(trend.100000034))

df102000002 <- df[which(df$site==102000002),]
max(df102000002$cDOY)-min(df102000002$cDOY)
g1<-gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df102000002)
trend.102000002 <- summary(gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df102000002))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df102000002 <- data.frame(site = 102000002, t(trend.102000002))

df104000058 <- df[which(df$site==104000058),]
max(df104000058$cDOY)-min(df104000058$cDOY)
g1<-gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df104000058)
trend.104000058 <- summary(gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df104000058))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df104000058 <- data.frame(site = 104000058, t(trend.104000058))

#order site results
tot_t<-rbind(trends, remies, zeries, zeries2, df100000034, df102000002, df104000058)
S_PIE_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("EPT_S_PIE_site", "EPT_S_PIE_Est", "EPT_S_PIE_SE", "EPT_S_PIE_p")
    colnames(S_PIE_df) <- xn
head(S_PIE_df)
nrow(S_PIE_df)

#########################################################################################################################
#ept rare.sppRich trend
hist(log10(ept$EPT_rare.sppRich))
#subset data for loop
df<-data.frame(ept$site_id_wMissing, ept$cday_of_year, ept$iYear, ept$cYear, log10(ept$EPT_rare.sppRich+1))
x <- c("site", "cDOY", "iYear", "cYear", "srr")
    colnames(df) <- x
head(df)
df2 <- df[!is.na(df$srr),]

#total values for each site to identify sites with no unique values over time
range <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  total.i = length(unique((sub$srr)))
  total.i <- data.frame(site = i, 
                        total.i)
  range <- rbind(range, total.i) ; rm(total.i, sub)
} ; rm(i)
head(range)
badi <- range[ which(range$total.i == 1),]
goodi <- range[ which(range$total.i > 1),]

# Create site subset for sites without all zeros
check = (df2$site %in% goodi$site) #make a list of T or F is there is a match in spp
df3 = df2[check,] #subset by this list

#calculate gls trends
trends <- NULL
for(i in unique(df3$site)){
  tryCatch({
    sub <- df3[df3$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(srr ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(srr ~ cDOY + cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[3, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
    }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)
head(trends)
##
 
Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn

df107000185 <- df[which(df$site==107000185),]
max(df107000185$cDOY)-min(df107000185$cDOY)
g1<-gls(srr ~ cYear + cDOY, na.action=na.omit, data = df107000185)
trend.107000185 <- summary(gls(srr ~ cYear + cDOY, na.action=na.omit, data = df107000185))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df107000185 <- data.frame(site = 107000185, t(trend.107000185))

df114000077 <- df[which(df$site==114000077),]
max(df114000077$cDOY)-min(df114000077$cDOY)
g1<-gls(srr ~ cYear + cDOY, na.action=na.omit, data = df114000077)
trend.114000077 <- summary(gls(srr ~ cYear + cDOY, na.action=na.omit, data = df114000077))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df114000077 <- data.frame(site = 114000077, t(trend.114000077))

#order site results
tot_t<-rbind(trends, zeries, df107000185, df114000077)
srr_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("EPT_rare.sppRich_site", "EPT_rare.sppRich_Est", "EPT_rare.sppRich_SE", "EPT_rare.sppRich_p")
    colnames(srr_df) <- xn
head(srr_df)
nrow(srr_df)

#########################################################################################################################

#########################################################################################################################
##############################################

glm.output<- cbind(sr_df,srr_df,sh_df,e10_df,ab_df,S_PIE_df,to_df)
head(glm.output)

glm.output$EPT_E10_site<-as.numeric(glm.output$EPT_E10_site)

identical(glm.output$EPT_SppRich_site,glm.output$EPT_rare.sppRich_site)
identical(glm.output$EPT_SppRich_site,glm.output$EPT_Shan_site)
identical(glm.output$EPT_SppRich_site,glm.output$EPT_E10_site)
identical(glm.output$EPT_SppRich_site,glm.output$EPT_Abund_site)
identical(glm.output$EPT_SppRich_site,glm.output$EPT_S_PIE_site)
identical(glm.output$EPT_SppRich_site,glm.output$EPT_Turnover_site)

library(dplyr)
new_df<-select(glm.output, -c("EPT_rare.sppRich_site","EPT_Shan_site","EPT_E10_site",
"EPT_Abund_site", "EPT_S_PIE_site", "EPT_Turnover_site"))

names(new_df)[names(new_df) == 'EPT_SppRich_site'] <- 'site'
head(new_df)

##
write.csv(new_df,"EPT_glmOutput.csv")
##

################################################################
###########################################################


