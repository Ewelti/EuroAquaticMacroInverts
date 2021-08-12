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
# 1 datasets with no insect in whole time series and should be removed for insect: 
#109000044
insect <- DATA2[which(DATA2$site_id!=109000044),]

attach(insect)
head(insect)
#####################################################################
#no sci notation
options(scipen = 999)
#################################################################
#insect spp_richness trend
hist(insect$insect_SppRich)

#subset data for loop
df<-data.frame(insect$site_id_wMissing, insect$cday_of_year, insect$iYear, insect$cYear, insect$insect_SppRich)
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

df103000121 <- df[which(df$site==103000121),]
max(df103000121$cDOY)-min(df103000121$cDOY)
g1<-gls(sr ~ cYear + cDOY, na.action=na.omit, data = df103000121)
trend.103000121 <- summary(gls(sr ~ cYear + cDOY, na.action=na.omit, data = df103000121))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df103000121 <- data.frame(site = 103000121, t(trend.103000121))

df109000411 <- df[which(df$site==109000411),]
max(df109000411$cDOY)-min(df109000411$cDOY)
g1<-gls(sr ~ cYear + cDOY, na.action=na.omit, data = df109000411)
trend.109000411 <- summary(gls(sr ~ cYear + cDOY, na.action=na.omit, data = df109000411))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000411 <- data.frame(site = 109000411, t(trend.109000411))

df111000061 <- df[which(df$site==111000061),]
max(df111000061$cDOY)-min(df111000061$cDOY)
g1<-gls(sr ~ cYear + cDOY, na.action=na.omit, data = df111000061)
trend.111000061 <- summary(gls(sr ~ cYear + cDOY, na.action=na.omit, data = df111000061))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #temporal autocorrelation :/
df111000061 <- data.frame(site = 111000061, Value = "NA", Std.Error = "NA", p.value = "NA") #NAs for trends from this site

#order site results
tot_t<-rbind(trends, df103000121, df109000411, df111000061)
sr_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("insect_SppRich_site", "insect_SppRich_Est", "insect_SppRich_SE", "insect_SppRich_p")
    colnames(sr_df) <- xn
head(sr_df)
nrow(sr_df)
################################################
################################################################################################
#insect abundance trend
hist(log10(insect$insect_Abund))

#subset data for loop
df<-data.frame(insect$site_id_wMissing, insect$cday_of_year, insect$iYear, insect$cYear, log10(insect$insect_Abund+1))
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

df104000122 <- df[which(df$site==104000122),]
max(df104000122$cDOY)-min(df104000122$cDOY)
g1<-gls(ab ~ cYear + cDOY, na.action=na.omit, data = df104000122)
trend.104000122 <- summary(gls(ab ~ cYear + cDOY, na.action=na.omit, data = df104000122))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df104000122 <- data.frame(site = 104000122, t(trend.104000122))

df107000085 <- df[which(df$site==107000085),]
max(df107000085$cDOY)-min(df107000085$cDOY)
g1<-gls(ab ~ cYear + cDOY, na.action=na.omit, data = df107000085)
trend.107000085 <- summary(gls(ab ~ cYear + cDOY, na.action=na.omit, data = df107000085))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #temporal autocorrelation
df107000085 <- data.frame(site = 107000085, Value = "NA", Std.Error = "NA", p.value = "NA") #NAs for trends from this site

#order site results
tot_t<-rbind(trends, df104000122, df107000085)
ab_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("insect_Abund_site", "insect_Abund_Est", "insect_Abund_SE", "insect_Abund_p")
    colnames(ab_df) <- xn
head(ab_df)
nrow(ab_df)

###################################################################################################################################
#insect shannon trend
hist(insect$insect_Shan)

#subset data for loop
df<-data.frame(insect$site_id_wMissing, insect$cday_of_year, insect$iYear, insect$cYear, insect$insect_Shan)
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

df121000096 <- df[which(df$site==121000096),]
max(df121000096$cDOY)-min(df121000096$cDOY)
g1<-gls(sh ~ cYear + cDOY, na.action=na.omit, data = df121000096)
trend.121000096 <- summary(gls(sh ~ cYear + cDOY, na.action=na.omit, data = df121000096))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df121000096 <- data.frame(site = 121000096, t(trend.121000096))

#order site results
tot_t<-rbind(trends, zeries, df121000096)
sh_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("insect_Shan_site", "insect_Shan_Est", "insect_Shan_SE", "insect_Shan_p")
    colnames(sh_df) <- xn
head(sh_df)
nrow(sh_df)
#########################################################################################################################
###################################################################################################################################
#insect E10 trend
hist(insect$insect_E10)
#subset data for loop
df<-data.frame(insect$site_id_wMissing, insect$cday_of_year, insect$iYear, insect$cYear, insect$insect_E10)
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

df107000204 <- df[which(df$site==107000204),]
max(df107000204$cDOY)-min(df107000204$cDOY)
g1<-gls(e10 ~ cYear + cDOY, na.action=na.omit, data = df107000204)
trend.107000204 <- summary(gls(e10 ~ cYear + cDOY, na.action=na.omit, data = df107000204))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df107000204 <- data.frame(site = 107000204, t(trend.107000204))

df108000046 <- df[which(df$site==108000046),]
max(df108000046$cDOY)-min(df108000046$cDOY)
g1<-gls(e10 ~ cYear + cDOY, na.action=na.omit, data = df108000046)
trend.108000046 <- summary(gls(e10 ~ cYear + cDOY, na.action=na.omit, data = df108000046))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000046 <- data.frame(site = 108000046, t(trend.108000046))

#order site results
tot_t<-rbind(trends, zeries, zeries2, df107000204, df108000046)
e10_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("insect_E10_site", "insect_E10_Est", "insect_E10_SE", "insect_E10_p")
    colnames(e10_df) <- xn
head(e10_df)
nrow(e10_df)
#########################################################################################################################
###################################################################################################################################
#insect Turnover trend
hist(insect$insect_Turnover)
#subset data for loop
df<-data.frame(insect$site_id_wMissing, insect$cday_of_year, insect$iYear, insect$cYear, insect$insect_Turnover)
x <- c("site", "cDOY", "iYear", "cYear", "Turnover")
    colnames(df) <- x
head(df)
df2 <- df[!is.na(df$Turnover),]

#total values for each site to identify sites with no unique values
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

df108000098 <- df[which(df$site==108000098),]
max(df108000098$cDOY)-min(df108000098$cDOY)
g1<-gls(Turnover ~ cYear + cDOY, na.action=na.omit, data = df108000098)
trend.108000098 <- summary(gls(Turnover ~ cYear+ cDOY, na.action=na.omit, data = df108000098))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000098 <- data.frame(site = 108000098, t(trend.108000098))

df108000128 <- df[which(df$site==108000128),]
max(df108000128$cDOY)-min(df108000128$cDOY)
g1<-gls(Turnover ~ cYear + cDOY, na.action=na.omit, data = df108000128)
trend.108000128 <- summary(gls(Turnover ~ cYear+ cDOY, na.action=na.omit, data = df108000128))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000128 <- data.frame(site = 108000128, t(trend.108000128))

df108000130 <- df[which(df$site==108000130),]
max(df108000130$cDOY)-min(df108000130$cDOY)
g1<-gls(Turnover ~ cYear + cDOY, na.action=na.omit, data = df108000130)
trend.108000130 <- summary(gls(Turnover ~ cYear+ cDOY, na.action=na.omit, data = df108000130))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000130 <- data.frame(site = 108000130, t(trend.108000130))

df109000067 <- df[which(df$site==109000067),]
max(df109000067$cDOY)-min(df109000067$cDOY)
g1<-gls(Turnover ~ cYear + cDOY, na.action=na.omit, data = df109000067)
trend.109000067 <- summary(gls(Turnover ~ cYear+ cDOY, na.action=na.omit, data = df109000067))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000067 <- data.frame(site = 109000067, t(trend.109000067))

df109000124 <- df[which(df$site==109000124),]
max(df109000124$cDOY)-min(df109000124$cDOY)
g1<-gls(Turnover ~ cYear + cDOY, na.action=na.omit, data = df109000124)
trend.109000124 <- summary(gls(Turnover ~ cYear+ cDOY, na.action=na.omit, data = df109000124))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000124 <- data.frame(site = 109000124, t(trend.109000124))

df109000318 <- df[which(df$site==109000318),]
max(df109000318$cDOY)-min(df109000318$cDOY)
g1<-gls(Turnover ~ cYear + cDOY, na.action=na.omit, data = df109000318)
trend.109000318 <- summary(gls(Turnover ~ cYear+ cDOY, na.action=na.omit, data = df109000318))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000318 <- data.frame(site = 109000318, t(trend.109000318))

#order site results
tot_t<-rbind(trends, zeries, zeries2, df108000098, df108000128, df108000130, df109000067, df109000124, df109000318)
to_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("insect_Turnover_site", "insect_Turnover_Est", "insect_Turnover_SE", "insect_Turnover_p")
    colnames(to_df) <- xn
head(to_df)
nrow(to_df)

###################################################################################################################################
#insect S_PIE trend
hist(log10(insect$insect_S_PIE))
#subset data for loop
df<-data.frame(insect$site_id_wMissing, insect$cday_of_year, insect$iYear, insect$cYear, log10(insect$insect_S_PIE+1))
x <- c("site", "cDOY", "iYear", "cYear", "S_PIE")
    colnames(df) <- x
head(df)
df2 <- df[!is.na(df$S_PIE),]

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
nrow(trends)
##
 
Value <- rep(0,each=nrow(badi))
bla <- rep("NA",each=nrow(badi))
zeries <- data.frame(badi$site, Value, bla, bla) 
xn <- c("site", "Value", "Std.Error", "p.value")
    colnames(zeries) <- xn
nrow(zeries)

bla <- rep("NA",each=nrow(badi2))
zeries2 <- data.frame(badi2$site, bla, bla, bla) 
colnames(zeries2) <- xn
nrow(zeries2)

df100000011 <- df[which(df$site==100000011),]
max(df100000011$cDOY)-min(df100000011$cDOY)
g1<-gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df100000011)
trend.100000011 <- summary(gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df100000011))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df100000011 <- data.frame(site = 100000011, t(trend.100000011))

df100000167 <- df[which(df$site==100000167),]
max(df100000167$cDOY)-min(df100000167$cDOY)
g1<-gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df100000167)
trend.100000167 <- summary(gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df100000167))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df100000167 <- data.frame(site = 100000167, t(trend.100000167))

df108000123 <- df[which(df$site==108000123),]
max(df108000123$cDOY)-min(df108000123$cDOY)
g1<-gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df108000123)
trend.108000123 <- summary(gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df108000123))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000123 <- data.frame(site = 108000123, t(trend.108000123))

df114000083 <- df[which(df$site==114000083),]
max(df114000083$cDOY)-min(df114000083$cDOY)
g1<-gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df114000083)
trend.114000083 <- summary(gls(S_PIE ~ cYear + cDOY, na.action=na.omit, data = df114000083))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df114000083 <- data.frame(site = 114000083, t(trend.114000083))

#order site results
tot_t<-rbind(trends, zeries, zeries2, df100000011, df100000167, df108000123, df114000083)
S_PIE_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("insect_S_PIE_site", "insect_S_PIE_Est", "insect_S_PIE_SE", "insect_S_PIE_p")
    colnames(S_PIE_df) <- xn
head(S_PIE_df)
nrow(S_PIE_df)

#########################################################################################################################
#insect rare.sppRich trend
hist(log10(insect$insect_rare.sppRich))
#subset data for loop
df<-data.frame(insect$site_id_wMissing, insect$cday_of_year, insect$iYear, insect$cYear, log10(insect$insect_rare.sppRich+1))
x <- c("site", "cDOY", "iYear", "cYear", "srr")
    colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
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

df103000590 <- df[which(df$site==103000590),]
max(df103000590$cDOY)-min(df103000590$cDOY)
g1<-gls(srr ~ cYear + cDOY, na.action=na.omit, data = df103000590)
trend.103000590 <- summary(gls(srr ~ cYear + cDOY, na.action=na.omit, data = df103000590))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df103000590 <- data.frame(site = 103000590, t(trend.103000590))

#order site results
tot_t<-rbind(trends, df103000590)
srr_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("insect_rare.sppRich_site", "insect_rare.sppRich_Est", "insect_rare.sppRich_SE", "insect_rare.sppRich_p")
    colnames(srr_df) <- xn
head(srr_df)
nrow(srr_df)

#########################################################################################################################

#########################################################################################################################
##############################################

glm.output<- cbind(sr_df,srr_df,sh_df,e10_df,ab_df,S_PIE_df,to_df)
head(glm.output)

glm.output$insect_E10_site<-as.numeric(glm.output$insect_E10_site)

identical(glm.output$insect_SppRich_site,glm.output$insect_rare.sppRich_site)
identical(glm.output$insect_SppRich_site,glm.output$insect_Shan_site)
identical(glm.output$insect_SppRich_site,glm.output$insect_E10_site)
identical(glm.output$insect_SppRich_site,glm.output$insect_Abund_site)
identical(glm.output$insect_SppRich_site,glm.output$insect_S_PIE_site)
identical(glm.output$insect_SppRich_site,glm.output$insect_Turnover_site)

library(dplyr)
new_df<-select(glm.output, -c("insect_rare.sppRich_site","insect_Shan_site","insect_E10_site",
"insect_Abund_site", "insect_S_PIE_site", "insect_Turnover_site"))

names(new_df)[names(new_df) == 'insect_SppRich_site'] <- 'site'
head(new_df)

##
write.csv(new_df,"insect_glmOutput.csv")
##

################################################################
###########################################################


