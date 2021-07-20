##Set working directory
#setwd("C:/Users/ewelti/Desktop/aquatic_data")

# load libraries
library(nlme)

# attach data
DATA1_list <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears.csv", header=T) # change file name according to the time series to be analyzed
DATA2 <- DATA1_list[!is.na(DATA1_list$site_id_wMissing),]
attach(DATA2)
head(DATA2)
turnover<-as.numeric(turnover) # make numeric

##############################make more date variables
#centre Year - helps model convergence to center variables for the model
DATA2$cYear <- DATA2$year_wMissing - median(DATA2$year_wMissing)

#or just have as an index starting from 1
DATA2$iYear <- DATA2$year_wMissing - min(DATA2$year_wMissing)+1

#centre day of year
DATA2$cday_of_year <- DATA2$day_of_year - median(DATA2$day_of_year,na.rm=T)

####################################
#no sci notation
options(scipen = 999)
#################################################################
#spp_richness trend
hist(spp_richness)

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, DATA2$spp_richness)
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

#order site results
sr_df <- trends[order(trends$site),] 

#rename gls output
xn <- c("SppRich_site", "SppRich_Est", "SppRich_SE", "SppRich_p")
colnames(sr_df) <- xn
head(sr_df)
nrow(sr_df)
mean(sr_df$SppRich_Est)
##########################################################################
#################################################################
#spp_richness rarefied trend
hist(spp_rich_rare)

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, DATA2$spp_rich_rare)
x <- c("site", "cDOY", "iYear", "cYear", "srra")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(srra ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(srra ~ cDOY + cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[3, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

df116000007 <- df[which(df$site==116000007),]
max(df116000007$cDOY)-min(df116000007$cDOY)
g1<-gls(srra ~ cYear + cDOY, na.action=na.omit, data = df116000007)
trend.116000007 <- summary(gls(srra ~ cYear + cDOY, na.action=na.omit, data = df116000007))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df116000007 <- data.frame(site = 116000007, t(trend.116000007))

tot_t<-rbind(trends, df116000007)
srra_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("SppRichRare_site", "SppRichRare_Est", "SppRichRare_SE", "SppRichRare_p")
colnames(srra_df) <- xn
head(srra_df)
nrow(srra_df)
mean(srra_df$SppRichRare_Est)
##########################################################################
#simpsonsD
hist(simpsonsD)

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, DATA2$simpsonsD)
x <- c("site", "cDOY", "iYear", "cYear", "si")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(si ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(si ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

#check problem sites
df117000040 <- df[which(df$site==117000040),]
max(df117000040$cDOY)-min(df117000040$cDOY)
g1<-gls(si ~ cYear + cDOY, na.action=na.omit, data = df117000040)
trend.117000040 <- summary(gls(si ~ cYear + cDOY, na.action=na.omit, data = df117000040))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df117000040 <- data.frame(site = 117000040, t(trend.117000040))

df121000223 <- df[which(df$site==121000223),]
max(df121000223$cDOY)-min(df121000223$cDOY)
g1<-gls(si ~ cYear + cDOY, na.action=na.omit, data = df121000223)
trend.121000223 <- summary(gls(si ~ cYear + cDOY, na.action=na.omit, data = df121000223))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df121000223 <- data.frame(site = 121000223, t(trend.121000223))

tot_t<-rbind(trends, df117000040, df121000223)
si_df <- tot_t[order(tot_t$site),] 

xn <- c("SimpD_site", "SimpD_Est", "SimpD_SE", "SimpD_p")
colnames(si_df) <- xn
head(si_df)
nrow(si_df)
mean(si_df$SimpD_Est)
##############################################
#shannonsH
hist(shannonsH)

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, DATA2$shannonsH)
x <- c("site", "cDOY", "iYear", "cYear", "sh")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(sh ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(sh ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

#check problem sites
df107000068 <- df[which(df$site==107000068),]
max(df107000068$cDOY)-min(df107000068$cDOY)
g1<-gls(sh ~ cYear + cDOY, na.action=na.omit, data = df107000068)
trend.107000068 <- summary(gls(sh ~ cYear + cDOY, na.action=na.omit, data = df107000068))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df107000068 <- data.frame(site = 107000068, t(trend.107000068))

tot_t<-rbind(trends, df107000068)
sh_df <- tot_t[order(tot_t$site),] 

xn <- c("ShanH_site", "ShanH_Est", "ShanH_SE", "ShanD_p")
colnames(sh_df) <- xn
head(sh_df)
nrow(sh_df)
mean(sh_df$ShanH_Est)
##############################################
#evennessJ
hist(evennessJ)

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, DATA2$evennessJ)
x <- c("site", "cDOY", "iYear", "cYear", "ev")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(ev ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(ev ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

#check problem sites
df100000075 <- df[which(df$site==100000075),]
max(df100000075$cDOY)-min(df100000075$cDOY)
g1<-gls(ev ~ cYear + cDOY, na.action=na.omit, data = df100000075)
trend.100000075 <- summary(gls(ev ~ cYear + cDOY, na.action=na.omit, data = df100000075))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df100000075 <- data.frame(site = 100000075, t(trend.100000075))

df108000084 <- df[which(df$site==108000084),]
max(df108000084$cDOY)-min(df108000084$cDOY)
g1<-gls(ev ~ cYear + cDOY, na.action=na.omit, data = df108000084)
trend.108000084 <- summary(gls(ev ~ cYear + cDOY, na.action=na.omit, data = df108000084))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000084 <- data.frame(site = 108000084, t(trend.108000084))

df121000223 <- df[which(df$site==121000223),]
max(df121000223$cDOY)-min(df121000223$cDOY)
g1<-gls(ev ~ cYear + cDOY, na.action=na.omit, data = df121000223)
trend.121000223 <- summary(gls(ev ~ cYear + cDOY, na.action=na.omit, data = df121000223))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df121000223 <- data.frame(site = 121000223, t(trend.121000223))

tot_t<-rbind(trends, df100000075, df108000084, df121000223)
ev_df <- tot_t[order(tot_t$site),] 

xn <- c("EvenJ_site", "EvenJ_Est", "EvenJ_SE", "EvenJ_p")
colnames(ev_df) <- xn
head(ev_df)
nrow(ev_df)
mean(ev_df$EvenJ_Est)
###############################################################
#abundance
hist(log10(abundance))

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, log10(DATA2$abundance))
x <- c("site", "cDOY", "iYear", "cYear", "abun")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(abun ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(abun ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

df109000279 <- df[which(df$site==109000279),]
max(df109000279$cDOY)-min(df109000279$cDOY)
g1<-gls(abun ~ cYear + cDOY, na.action=na.omit, data = df109000279)
trend.109000279 <- summary(gls(abun ~ cYear + cDOY, na.action=na.omit, data = df109000279))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000279 <- data.frame(site = 109000279, t(trend.109000279))

tot_t<-rbind(trends, df109000279)
ab_df <- tot_t[order(tot_t$site),] 

xn <- c("Abun_site", "Abun_Est", "Abun_SE", "Abun_p")
colnames(ab_df) <- xn
nrow(ab_df)

mean(ab_df$Abun_Est)
##############################################
#################################################################
#S_PIE trend
hist(log10(S_PIE))

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, log10(DATA2$S_PIE))
x <- c("site", "cDOY", "iYear", "cYear", "spie")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(spie ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(spie ~ cDOY + cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[3, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

df108000027 <- df[which(df$site==108000027),]
max(df108000027$cDOY)-min(df108000027$cDOY)
g1<-gls(spie ~ cYear + cDOY, na.action=na.omit, data = df108000027)
trend.108000027 <- summary(gls(spie ~ cYear + cDOY, na.action=na.omit, data = df108000027))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000027 <- data.frame(site = 108000027, t(trend.108000027))

tot_t<-rbind(trends, df108000027)
spie_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("S_PIE_site", "S_PIE_Est", "S_PIE_SE", "S_PIE_p")
colnames(spie_df) <- xn
head(spie_df)
nrow(spie_df)
mean(spie_df$S_PIE_Est)
#####################################################################################
#turnover
hist(turnover)

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, DATA2$turnover)
x <- c("site", "cDOY", "iYear", "cYear", "to")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(to ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(to ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

#check temporal autocorrelation in problem sites
df108000087 <- df[which(df$site==108000087),]
max(df108000087$cDOY)-min(df108000087$cDOY)
g1<-gls(to ~ cYear + cDOY, na.action=na.omit, data = df108000087)
trend.108000087 <- summary(gls(to ~ cYear + cDOY, na.action=na.omit, data = df108000087))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000087 <- data.frame(site = 108000087, t(trend.108000087))

df108000131 <- df[which(df$site==108000131),]
max(df108000131$cDOY)-min(df108000131$cDOY)
g1<-gls(to ~ cYear + cDOY, na.action=na.omit, data = df108000131)
trend.108000131 <- summary(gls(to ~ cYear + cDOY, na.action=na.omit, data = df108000131))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000131 <- data.frame(site = 108000131, t(trend.108000131))

df109000128 <- df[which(df$site==109000128),]
max(df109000128$cDOY)-min(df109000128$cDOY)
g1<-gls(to ~ cYear + cDOY, na.action=na.omit, data = df109000128)
trend.109000128 <- summary(gls(to ~ cYear + cDOY, na.action=na.omit, data = df109000128))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000128 <- data.frame(site = 109000128, t(trend.109000128))

df116000004 <- df[which(df$site==116000004),]
max(df116000004$cDOY)-min(df116000004$cDOY)
g1<-gls(to ~ cYear + cDOY, na.action=na.omit, data = df116000004)
trend.116000004 <- summary(gls(to ~ cYear + cDOY, na.action=na.omit, data = df116000004))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df116000004 <- data.frame(site = 116000004, t(trend.116000004))

df121000013 <- df[which(df$site==121000013),]
max(df121000013$cDOY)-min(df121000013$cDOY)
g1<-gls(to ~ cYear + cDOY, na.action=na.omit, data = df121000013)
trend.121000013 <- summary(gls(to ~ cYear + cDOY, na.action=na.omit, data = df121000013))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df121000013 <- data.frame(site = 121000013, t(trend.121000013))

df121000193 <- df[which(df$site==121000193),]
max(df121000193$cDOY)-min(df121000193$cDOY)
g1<-gls(to ~ cYear + cDOY, na.action=na.omit, data = df121000193)
trend.121000193 <- summary(gls(to ~ cYear + cDOY, na.action=na.omit, data = df121000193))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df121000193 <- data.frame(site = 121000193, t(trend.121000193))

#bind trend results
tot_t<-rbind(trends, df108000087, df108000131, df109000128, df116000004, df121000013, df121000193)
to_df <- tot_t[order(tot_t$site),] 

xn <- c("TurnO_site", "TurnO_Est", "TurnO_SE", "TurnO_p")
colnames(to_df) <- xn
nrow(to_df)

mean(to_df$TurnO_Est)
##############################################
#F_to
hist(F_to)

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, DATA2$F_to)
x <- c("site", "cDOY", "iYear", "cYear", "F_to")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(F_to ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(F_to ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

df103000552 <- df[which(df$site==103000552),]
df103000552$F_to
df103000552 <- data.frame(site = 103000552, Value = 0, Std.Error = "NA", p.value = "NA")

df104000142 <- df[which(df$site==104000142),]
df104000142$F_to
df104000142 <- data.frame(site = 104000142, Value = 0, Std.Error = "NA", p.value = "NA")

tot_t<-rbind(trends, df103000552, df104000142)
F_to_df <- tot_t[order(tot_t$site),] 

xn <- c("F_to_site", "F_to_Est", "F_to_SE", "F_to_p")
colnames(F_to_df) <- xn
nrow(F_to_df)

mean(F_to_df$F_to_Est)
##############################################
#FRic
hist(log10(FRic))

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, log10(DATA2$FRic))
x <- c("site", "cDOY", "iYear", "cYear", "FR")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(FR ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(FR ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

df103000559 <- df[which(df$site==103000559),]
max(df103000559$cDOY)-min(df103000559$cDOY)
g1<-gls(FR ~ cYear, na.action=na.omit, data = df103000559)
trend.103000559 <- summary(gls(FR ~ cYear, na.action=na.omit, data = df103000559))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df103000559 <- data.frame(site = 103000559, t(trend.103000559))

df109000348 <- df[which(df$site==109000348),]
max(df109000348$cDOY)-min(df109000348$cDOY)
g1<-gls(FR ~ cYear + cDOY, na.action=na.omit, data = df109000348)
trend.109000348 <- summary(gls(FR ~ cYear + cDOY, na.action=na.omit, data = df109000348))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000348 <- data.frame(site = 109000348, t(trend.109000348))

tot_t<-rbind(trends, df103000559, df109000348)
FRic_df <- tot_t[order(tot_t$site),] 

xn <- c("FRic_site", "FRic_Est", "FRic_SE", "FRic_p")
colnames(FRic_df) <- xn
nrow(FRic_df)

mean(FRic_df$FRic_Est)
##############################################
#FEve
hist(FEve)

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, DATA2$FEve)
x <- c("site", "cDOY", "iYear", "cYear", "FEve")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(FEve ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(FEve ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

df100000297 <- df[which(df$site==100000297),]
max(df100000297$cDOY)-min(df100000297$cDOY)
g1<-gls(FEve ~ cYear + cDOY, na.action=na.omit, data = df100000297)
trend.100000297 <- summary(gls(FEve ~ cYear + cDOY, na.action=na.omit, data = df100000297))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df100000297 <- data.frame(site = 100000297, t(trend.100000297))

df109000006 <- df[which(df$site==109000006),]
max(df109000006$cDOY)-min(df109000006$cDOY)
g1<-gls(FEve ~ cYear + cDOY, na.action=na.omit, data = df109000006)
trend.109000006 <- summary(gls(FEve ~ cYear + cDOY, na.action=na.omit, data = df109000006))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000006 <- data.frame(site = 109000006, t(trend.109000006))

df121000004 <- df[which(df$site==121000004),]
max(df121000004$cDOY)-min(df121000004$cDOY)
g1<-gls(FEve ~ cYear + cDOY, na.action=na.omit, data = df121000004)
trend.121000004 <- summary(gls(FEve ~ cYear + cDOY, na.action=na.omit, data = df121000004))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df121000004 <- data.frame(site = 121000004, t(trend.121000004))

tot_t<-rbind(trends, df100000297, df109000006, df121000004)
FEve_df <- tot_t[order(tot_t$site),] 

xn <- c("FEve_site", "FEve_Est", "FEve_SE", "FEve_p")
colnames(FEve_df) <- xn
nrow(FEve_df)

mean(FEve_df$FEve_Est)
##############################################
#FDiv
hist(FDiv)

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, DATA2$FDiv)
x <- c("site", "cDOY", "iYear", "cYear", "FDiv")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(FDiv ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(FDiv ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

df103000590 <- df[which(df$site==103000590),]
max(df103000590$cDOY)-min(df103000590$cDOY)
g1<-gls(FDiv ~ cYear + cDOY, na.action=na.omit, data = df103000590)
trend.103000590 <- summary(gls(FDiv ~ cYear + cDOY, na.action=na.omit, data = df103000590))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df103000590 <- data.frame(site = 103000590, t(trend.103000590))

df108000059 <- df[which(df$site==108000059),]
max(df108000059$cDOY)-min(df108000059$cDOY)
g1<-gls(FDiv ~ cYear + cDOY, na.action=na.omit, data = df108000059)
trend.108000059 <- summary(gls(FDiv ~ cYear + cDOY, na.action=na.omit, data = df108000059))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000059 <- data.frame(site = 108000059, t(trend.108000059))

tot_t<-rbind(trends, df103000590, df108000059)
FDiv_df <- tot_t[order(tot_t$site),] 

xn <- c("FDiv_site", "FDiv_Est", "FDiv_SE", "FDiv_p")
colnames(FDiv_df) <- xn
nrow(FDiv_df)

mean(FDiv_df$FDiv_Est)
##############################################
#FDis
hist(FDis)

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, DATA2$FDis)
x <- c("site", "cDOY", "iYear", "cYear", "FDis")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(FDis ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(FDis ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

df108000027 <- df[which(df$site==108000027),]
max(df108000027$cDOY)-min(df108000027$cDOY)
g1<-gls(FDis ~ cYear + cDOY, na.action=na.omit, data = df108000027)
trend.108000027 <- summary(gls(FDis ~ cYear + cDOY, na.action=na.omit, data = df108000027))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000027 <- data.frame(site = 108000027, t(trend.108000027))

df109000367 <- df[which(df$site==109000367),]
max(df109000367$cDOY)-min(df109000367$cDOY)
g1<-gls(FDis ~ cYear + cDOY, na.action=na.omit, data = df109000367)
trend.109000367 <- summary(gls(FDis ~ cYear + cDOY, na.action=na.omit, data = df109000367))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000367 <- data.frame(site = 109000367, t(trend.109000367))

tot_t<-rbind(trends, df108000027, df109000367)
FDis_df <- tot_t[order(tot_t$site),] 

xn <- c("FDis_site", "FDis_Est", "FDis_SE", "FDis_p")
colnames(FDis_df) <- xn
nrow(FDis_df)

mean(FDis_df$FDis_Est)
##############################################
#RaoQ
hist(RaoQ)

#subset data for loop
df<-data.frame(DATA2$site_id_wMissing, DATA2$cday_of_year, DATA2$iYear, DATA2$cYear, DATA2$RaoQ)
x <- c("site", "cDOY", "iYear", "cYear", "RaoQ")
colnames(df) <- x
head(df)

#calculate gls trends
trends <- NULL
for(i in unique(df$site)){
  tryCatch({
    sub <- df[df$site == i, ]
    maxDiffDays = max(sub$cDOY)-min(sub$cDOY)
    if(maxDiffDays < 30) {
      trend.i <- summary(gls(RaoQ ~ cYear, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    } else{
      trend.i <- summary(gls(RaoQ ~ cYear + cDOY, correlation = corAR1(form = ~ iYear),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    }
    trend.i <- data.frame(site = i, 
                          t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub, maxDiffDays)
  }, error=function(e){cat(unique(sub$site),conditionMessage(e), "\n")})
} ; rm(i)

df109000050 <- df[which(df$site==109000050),]
max(df109000050$cDOY)-min(df109000050$cDOY)
g1<-gls(RaoQ ~ cYear + cDOY, na.action=na.omit, data = df109000050)
trend.109000050 <- summary(gls(RaoQ ~ cYear + cDOY, na.action=na.omit, data = df109000050))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000050 <- data.frame(site = 109000050, t(trend.109000050))

tot_t<-rbind(trends, df109000050)
RaoQ_df <- tot_t[order(tot_t$site),]

xn <- c("RaoQ_site", "RaoQ_Est", "RaoQ_SE", "RaoQ_p")
colnames(RaoQ_df) <- xn
nrow(RaoQ_df)
mean(RaoQ_df$RaoQ_Est)
##############################################

glm.output<- cbind(sr_df,srra_df,si_df,sh_df,ev_df,ab_df,spie_df,to_df,F_to_df,FRic_df,FEve_df,FDiv_df,RaoQ_df)
head(glm.output)

glm.output$Abun_site<-as.numeric(glm.output$Abun_site)
glm.output$ShanH_site<-as.numeric(glm.output$ShanH_site)
glm.output$EvenJ_site<-as.numeric(glm.output$EvenJ_site)
glm.output$FDiv_site<-as.numeric(glm.output$FDiv_site)
glm.output$RaoQ_site<-as.numeric(glm.output$RaoQ_site)
glm.output$SppRich_site<-as.numeric(glm.output$SppRich_site)

identical(glm.output$SppRich_site,glm.output$SppRichRare_site)
identical(glm.output$SppRich_site,glm.output$SimpD_site)
identical(glm.output$SppRich_site,glm.output$ShanH_site)
identical(glm.output$SppRich_site,glm.output$EvenJ_site)
identical(glm.output$SppRich_site,glm.output$Abun_site)
identical(glm.output$SppRich_site,glm.output$S_PIE_site)
identical(glm.output$SppRich_site,glm.output$TurnO_site)
identical(glm.output$SppRich_site,glm.output$F_to_site)
identical(glm.output$SppRich_site,glm.output$FRic_site)
identical(glm.output$SppRich_site,glm.output$FEve_site)
identical(glm.output$SppRich_site,glm.output$FDiv_site)
identical(glm.output$SppRich_site,glm.output$RaoQ_site)

library(dplyr)
new_df<-select(glm.output, -c("SppRichRare_site","SimpD_site","ShanH_site","EvenJ_site","Abun_site", "S_PIE_site",
                              "TurnO_site","F_to_site","FRic_site","FEve_site","FDiv_site","RaoQ_site"))

names(new_df)[names(new_df) == 'SppRich_site'] <- 'site'
head(new_df)

##
write.csv(new_df,"glmOutput.csv")
##

################################################################
###########################################################
###################################################
##############################################