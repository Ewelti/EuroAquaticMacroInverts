##Set working directory
#setwd("C:/Users/ewelti/Desktop/aquatic_data")

###################################################
##############################################
# load libraries
library(nlme)

#####ALIENS
# attach data
#this data is the same as raw-data/All_indices_benthicMacroInverts_AllYears.csv
#accept that the sites IDed to family level are removed because non-native species are not able to be be determined at family level
#and sites with all no non-native species are removed because these make singular models
DATA3 <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienSub.csv", header=T) # change file name according to the time series to be analyzed
DATA4 <- DATA3[!is.na(DATA3$site_id_wMissing),]
attach(DATA4)
head(DATA4)
####################################
#no sci notation
options(scipen = 999)
#####################################
##abundProportion_alienSpp

#subset data for loop
df<-data.frame(DATA4$site_id_wMissing, DATA4$year_wMissing, DATA4$abundProportion_alienSpp)
x <- c("site", "yr", "PropAlienAbun")
    colnames(df) <- x
head(df)

unique(df$site)

df2 <- df[which(df$site!=107000164& df$site!=109000041& df$site!=114000091& df$site!=121000190),]

trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(PropAlienAbun ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

df107000164 <- df[which(df$site==107000164),]
g1<-gls(PropAlienAbun ~ yr,na.action=na.omit, data = df107000164)
trend.107000164 <- summary(gls(PropAlienAbun ~ yr,na.action=na.omit, data = df107000164))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #temporal autocorrelation :/
df107000164 <- data.frame(site = 107000164, Value = "NA", Std.Error = "NA", p.value = "NA")

df109000041 <- df[which(df$site==109000041),]
g1<-gls(PropAlienAbun ~ yr,na.action=na.omit, data = df109000041)
trend.109000041 <- summary(gls(PropAlienAbun ~ yr,na.action=na.omit, data = df109000041))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df109000041 <- data.frame(site = 109000041, t(trend.109000041))

df114000091 <- df[which(df$site==109000041),]
g1<-gls(PropAlienAbun ~ yr,na.action=na.omit, data = df114000091)
trend.114000091 <- summary(gls(PropAlienAbun ~ yr,na.action=na.omit, data = df114000091))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #temporal autocorrelation :/
df114000091 <- data.frame(site = 114000091, t(trend.114000091))

df121000190 <- df[which(df$site==121000190),]
g1<-gls(PropAlienAbun ~ yr,na.action=na.omit, data = df121000190)
trend.121000190 <- summary(gls(PropAlienAbun ~ yr,na.action=na.omit, data = df121000190))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df121000190 <- data.frame(site = 121000190, t(trend.121000190))

tot_t<-rbind(trends, df107000164, df109000041, df114000091, df121000190)

AlienAbunProp_df <- tot_t[order(tot_t$site),] 
xn <- c("AlienAbunProp_site", "AlienAbunProp_Est", "AlienAbunProp_SE", "AlienAbunProp_p")
    colnames(AlienAbunProp_df) <- xn
nrow(AlienAbunProp_df)

##############################################

##SppRichProportion_alienSpp
#subset data for loop
df<-data.frame(DATA4$site_id_wMissing, DATA4$year_wMissing, DATA4$SppRichProportion_alienSpp)
x <- c("site", "yr", "PropAlienSppRich")
    colnames(df) <- x
head(df)

unique(df$site)

trends <- NULL
for(i in unique(df$site)){
  sub <- df[df$site == i, ]
  trend.i <- summary(gls(PropAlienSppRich ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

AlienSppRichProp_df <- trends[order(trends$site),] 
xn <- c("AlienSppRichProp_site", "AlienSppRichProp_Est", "AlienSppRichProp_SE", "AlienSppRichProp_p")
    colnames(AlienSppRichProp_df) <- xn
nrow(AlienSppRichProp_df)

##############################################
##############################################

##abund_alienSpp
#subset data for loop
hist(DATA4$abund_alienSpp)

df<-data.frame(DATA4$site_id_wMissing, DATA4$year_wMissing, log10(DATA4$abund_alienSpp+1))
x <- c("site", "yr", "AlienAbun")
    colnames(df) <- x
head(df)

unique(df$site)

#df2 <- df[which(df$site!=108000057& df$site!=108000149& df$site!=109000401),] #withoutlog
df2 <- df[which(df$site!=107000057& df$site!=107000234& df$site!=109000169),]

trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(AlienAbun ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

df107000057 <- df[which(df$site==107000057),]
g1<-gls(AlienAbun ~ yr,na.action=na.omit, data = df107000057)
trend.107000057 <- summary(gls(AlienAbun ~ yr,na.action=na.omit, data = df107000057))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df107000057 <- data.frame(site = 107000057, t(trend.107000057))

df107000234 <- df[which(df$site==107000234),]
g1<-gls(AlienAbun ~ yr,na.action=na.omit, data = df107000234)
trend.107000234 <- summary(gls(AlienAbun ~ yr,na.action=na.omit, data = df107000234))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df107000234 <- data.frame(site = 107000234, t(trend.107000234))

df109000169 <- df[which(df$site==109000169),]
g1<-gls(AlienAbun ~ yr,na.action=na.omit, data = df109000169)
trend.109000169 <- summary(gls(AlienAbun ~ yr,na.action=na.omit, data = df109000169))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000169 <- data.frame(site = 109000169, t(trend.109000169))

tot_t<-rbind(trends, df107000057, df107000234, df109000169)

AlienAbun_df <- tot_t[order(tot_t$site),] 
xn <- c("AlienAbun_site", "AlienAbun_Est", "AlienAbun_SE", "AlienAbun_p")
    colnames(AlienAbun_df) <- xn
nrow(AlienAbun_df)

##############################################
##############################################

##SppRich_alienSpp
#subset data for loop
hist(DATA4$SppRich_alienSpp)

df<-data.frame(DATA4$site_id_wMissing, DATA4$year_wMissing, DATA4$SppRich_alienSpp)
x <- c("site", "yr", "AlienSppRich")
    colnames(df) <- x
head(df)

unique(DATA4$site_id_wMissing)

df2 <- df[which(df$site!=105000006& df$site!=105000012& df$site!=107000057& df$site!=107000126& 
df$site!=107000228& df$site!=108000057& df$site!=109000079& df$site!=109000080& df$site!=109000143&
df$site!=109000163& df$site!=109000190& df$site!=109000343& df$site!=109000374& df$site!=109000375&
df$site!=114000090),]

trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(AlienSppRich ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

df105000006 <- df[which(df$site==105000006),]
df105000006
df105000006 <- data.frame(site = 105000006, Value = "0", Std.Error = "0", p.value = "NA")

df105000012 <- df[which(df$site==105000012),]
df105000012
df105000012 <- data.frame(site = 105000012, Value = "0", Std.Error = "0", p.value = "NA")

df107000057 <- df[which(df$site==107000057),]
g1<-gls(AlienSppRich ~ yr,na.action=na.omit, data = df107000057)
trend.107000057 <- summary(gls(AlienSppRich ~ yr,na.action=na.omit, data = df107000057))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df107000057 <- data.frame(site = 107000057, t(trend.107000057))

df107000126 <- df[which(df$site==107000126),]
g1<-gls(AlienSppRich ~ yr,na.action=na.omit, data = df107000126)
trend.107000126 <- summary(gls(AlienSppRich ~ yr,na.action=na.omit, data = df107000126))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #temporal autocorrelation :/
df107000126 <- data.frame(site = 107000126, Value = "NA", Std.Error = "NA", p.value = "NA")

df107000228 <- df[which(df$site==107000228),]
g1<-gls(AlienSppRich ~ yr,na.action=na.omit, data = df107000228)
trend.107000228 <- summary(gls(AlienSppRich ~ yr,na.action=na.omit, data = df107000228))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df107000228 <- data.frame(site = 107000228, t(trend.107000228))

df108000057 <- df[which(df$site==108000057),]
g1<-gls(AlienSppRich ~ yr,na.action=na.omit, data = df108000057)
trend.108000057 <- summary(gls(AlienSppRich ~ yr,na.action=na.omit, data = df108000057))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000057 <- data.frame(site = 108000057, t(trend.108000057))

df109000079 <- df[which(df$site==109000079),]
g1<-gls(AlienSppRich ~ yr,na.action=na.omit, data = df109000079)
trend.109000079 <- summary(gls(AlienSppRich ~ yr,na.action=na.omit, data = df109000079))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000079 <- data.frame(site = 109000079, t(trend.109000079))

df109000080 <- df[which(df$site==109000080),]
df109000080
df109000080 <- data.frame(site = 109000080, Value = "0", Std.Error = "0", p.value = "NA")

df109000143 <- df[which(df$site==109000143),]
df109000143
df109000143 <- data.frame(site = 109000143, Value = "0", Std.Error = "0", p.value = "NA")

df109000163 <- df[which(df$site==109000163),]
df109000163
df109000163 <- data.frame(site = 109000163, Value = "0", Std.Error = "0", p.value = "NA")

df109000190 <- df[which(df$site==109000190),]
df109000190
df109000190 <- data.frame(site = 109000190, Value = "0", Std.Error = "0", p.value = "NA")

df109000343 <- df[which(df$site==109000343),]
df109000343
df109000343 <- data.frame(site = 109000343, Value = "0", Std.Error = "0", p.value = "NA")

df109000374 <- df[which(df$site==109000374),]
df109000374
df109000374 <- data.frame(site = 109000374, Value = "0", Std.Error = "0", p.value = "NA")

df109000375 <- df[which(df$site==109000375),]
df109000375
df109000375 <- data.frame(site = 109000375, Value = "0", Std.Error = "0", p.value = "NA")

df114000090 <- df[which(df$site==114000090),]
df114000090
df114000090 <- data.frame(site = 114000090, Value = "0", Std.Error = "0", p.value = "NA")

tot_t<-rbind(trends, df105000006, df105000012, df107000057, df107000126, df107000228, 
df108000057, df109000079, df109000080, df109000143, df109000163, 
df109000190, df109000343, df109000374, df109000375, df114000090)

AlienSppRich_df <- tot_t[order(tot_t$site),] 
xn <- c("AlienSppRich_site", "AlienSppRich_Est", "AlienSppRich_SE", "AlienSppRich_p")
    colnames(AlienSppRich_df) <- xn
nrow(AlienSppRich_df)

##############################################
##############################################

##abund_nativeSpp
#subset data for loop
hist(DATA4$abund_nativeSpp)

df<-data.frame(DATA4$site_id_wMissing, DATA4$year_wMissing, log10(DATA4$abund_nativeSpp+1))
x <- c("site", "yr", "nativeAbun")
    colnames(df) <- x
head(df)

unique(DATA4$site_id_wMissing)

df2 <- df[which(df$site!=109000037& df$site!=117000041),]

trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(nativeAbun ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

df109000037 <- df[which(df$site==109000037),]
g1<-gls(nativeAbun ~ yr,na.action=na.omit, data = df109000037)
trend.109000037 <- summary(gls(nativeAbun ~ yr,na.action=na.omit, data = df109000037))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #temporal autocorrelation :/
df109000037 <- data.frame(site = 109000037, Value = "NA", Std.Error = "NA", p.value = "NA")

df117000041 <- df[which(df$site==117000041),]
g1<-gls(nativeAbun ~ yr,na.action=na.omit, data = df117000041)
trend.117000041 <- summary(gls(nativeAbun ~ yr,na.action=na.omit, data = df117000041))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df117000041 <- data.frame(site = 117000041, t(trend.117000041))

tot_t<-rbind(trends, df109000037, df117000041)

nativeAbun_df <- tot_t[order(tot_t$site),] 
xn <- c("nativeAbun_site", "nativeAbun_Est", "nativeAbun_SE", "nativeAbun_p")
    colnames(nativeAbun_df) <- xn
nrow(nativeAbun_df)

##############################################
##############################################

##SppRich_nativeSpp
#subset data for loop
hist(DATA4$SppRich_nativeSpp)

df<-data.frame(DATA4$site_id_wMissing, DATA4$year_wMissing, DATA4$SppRich_nativeSpp)
x <- c("site", "yr", "nativeSppRich")
    colnames(df) <- x
head(df)

trends <- NULL
for(i in unique(df$site)){
  sub <- df[df$site == i, ]
  trend.i <- summary(gls(nativeSppRich ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

nativeSppRich_df <- trends[order(trends$site),] 
xn <- c("nativeSppRich_site", "nativeSppRich_Est", "nativeSppRich_SE", "nativeSppRich_p")
    colnames(nativeSppRich_df) <- xn
nrow(nativeSppRich_df)

##
##############################################

glm.output<- cbind(AlienAbunProp_df,AlienSppRichProp_df,AlienAbun_df,AlienSppRich_df,nativeAbun_df,nativeSppRich_df)
head(glm.output)

glm.output$AlienAbunProp_site<-as.numeric(glm.output$AlienAbunProp_site)
glm.output$AlienSppRichProp_site<-as.numeric(glm.output$AlienSppRichProp_site)
glm.output$AlienAbun_site<-as.numeric(glm.output$AlienAbun_site)
glm.output$AlienAlienSppRich_site<-as.numeric(glm.output$AlienAlienSppRich_site)
glm.output$nativeAbun_site<-as.numeric(glm.output$nativeAbun_site)
glm.output$nativeSppRich_site<-as.numeric(glm.output$nativeSppRich_site)

identical(glm.output$AlienAbunProp_site,glm.output$AlienSppRichProp_site)
identical(glm.output$AlienAbunProp_site,glm.output$AlienAbun_site)
identical(glm.output$AlienAbunProp_site,glm.output$AlienSppRich_site)
identical(glm.output$AlienAbunProp_site,glm.output$nativeAbun_site)
identical(glm.output$AlienAbunProp_site,glm.output$nativeSppRich_site)

library(dplyr)
new_df<-select(glm.output, -c("AlienSppRichProp_site","AlienAbun_site","AlienSppRich_site",
"nativeAbun_site","nativeSppRich_site"))

names(new_df)[names(new_df) == 'AlienAbunProp_site'] <- 'site'
head(new_df)

##
write.csv(new_df,"glmOut_AlienSlopes.csv")
##
