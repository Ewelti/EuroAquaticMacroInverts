##Set working directory
setwd("C:/Users/ewelti/Desktop/aquatic_data")

# load libraries
library(nlme)

# attach data
DATA1_list <- read.csv("All_indices_benthicMacroInverts.csv", header=T) # change file name according to the time series to be analyzed
attach(DATA1_list)
head(DATA1_list)
turnover<-as.numeric(turnover) # make numeric

####################################
#no sci notation
options(scipen = 999)
#################################################################
#spp_richness trend
hist((spp_richness))

#subset data for loop
df<-data.frame(site_id, year.1,spp_richness)
x <- c("site", "yr", "sr")
    colnames(df) <- x
head(df)

#remove problem site
df2 <- df[which(df$site!=108000041),] #df2 <- df[which(df$site!=100000293),]# if log transformed

#calculate gls trends
trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(sr ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

#check temporal autocorrelation in problem site
df108000041 <- df[which(df$site==108000041),]
g1<-gls(sr ~ yr,na.action=na.omit, data = df108000041)
trend.108000041 <- summary(gls(sr ~ yr,na.action=na.omit, data = df108000041))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #temporal autocorrelation :/
df108000041 <- data.frame(site = 108000041, Value = "NA", Std.Error = "NA", p.value = "NA") #NAs for trends from this site

#bind site results
tot_t<-rbind(trends, df108000041)
sr_df <- tot_t[order(tot_t$site),] 

#rename gls output
xn <- c("SppRich_site", "SppRich_Est", "SppRich_SE", "SppRich_p")
    colnames(sr_df) <- xn
head(sr_df)
##########################################################################
#simpsonsD
hist(log10(simpsonsD))
hist((simpsonsD))

#subset data for loop
df<-data.frame(site_id, year.1, simpsonsD)
x <- c("site", "yr", "si")
    colnames(df) <- x
head(df)

df2 <- df[which(df$site!=100000270& df$site!=103000631&
df$site!=103000632 & df$site!=105000004& df$site!=105000005& 
df$site!=109000011& df$site!=121000189),]

trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(si ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)


df100000270 <- df[which(df$site==100000270),]
g1<-gls(si ~ yr,na.action=na.omit, data = df100000270)
trend.100000270 <- summary(gls(si ~ yr,na.action=na.omit, data = df100000270))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df100000270 <- data.frame(site = 100000270, t(trend.100000270))

df103000631 <- df[which(df$site==103000631),]
g1<-gls(si ~ yr,na.action=na.omit, data = df103000631)
trend.103000631 <- summary(gls(si ~ yr,na.action=na.omit, data = df103000631))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #temporal autocorrelation :/
df103000631 <- data.frame(site = 103000631, Value = "NA", Std.Error = "NA", p.value = "NA")

df103000632 <- df[which(df$site==103000632),]
g1<-gls(si ~ yr,na.action=na.omit, data = df103000632)
trend.103000632 <- summary(gls(si ~ yr,na.action=na.omit, data = df103000632))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df103000632  <- data.frame(site = 103000632, t(trend.103000632))

df105000004 <- df[which(df$site==105000004),]
g1<-gls(si ~ yr,na.action=na.omit, data = df105000004)
trend.105000004 <- summary(gls(si ~ yr,na.action=na.omit, data = df105000004))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df105000004 <- data.frame(site = 105000004, t(trend.105000004))

df105000005 <- df[which(df$site==105000005),]
g1<-gls(si ~ yr,na.action=na.omit, data = df105000005)
trend.105000005 <- summary(gls(si ~ yr,na.action=na.omit, data = df105000005))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df105000005 <- data.frame(site = 105000005, t(trend.105000005))

df109000011 <- df[which(df$site==109000011),]
g1<-gls(si ~ yr,na.action=na.omit, data = df109000011)
trend.109000011 <- summary(gls(si ~ yr,na.action=na.omit, data = df109000011))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000011  <- data.frame(site = 109000011, t(trend.109000011))

df121000189 <- df[which(df$site==121000189),]
g1<-gls(si ~ yr,na.action=na.omit, data = df121000189)
trend.121000189 <- summary(gls(si ~ yr,na.action=na.omit, data = df121000189))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df121000189  <- data.frame(site = 121000189, t(trend.121000189))

tot_t<-rbind(trends, df100000270,df103000631, df103000632, df105000004, df105000005, df109000011, df121000189)
si_df <- tot_t[order(tot_t$site),] 

xn <- c("SimpD_site", "SimpD_Est", "SimpD_SE", "SimpD_p")
    colnames(si_df) <- xn
head(si_df)

##############################################
#shannonsH
hist((shannonsH))

#subset data for loop
df<-data.frame(site_id, year.1, shannonsH)
x <- c("site", "yr", "sh")
    colnames(df) <- x
head(df)

trends <- NULL
for(i in unique(df$site)){
  sub <- df[df$site == i, ]
  trend.i <- summary(gls(sh ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

sh_df <- trends[order(trends$site),] 
xn <- c("ShanH_site", "ShanH_Est", "ShanH_SE", "ShanD_p")
    colnames(sh_df) <- xn
head(sh_df)

##############################################
#evennessJ
hist((evennessJ))

#subset data for loop
df<-data.frame(site_id, year.1, evennessJ)
x <- c("site", "yr", "ev")
    colnames(df) <- x
head(df)

trends <- NULL
for(i in unique(df$site)){
  sub <- df[df$site == i, ]
  trend.i <- summary(gls(ev ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

ev_df <- trends[order(trends$site),] 
xn <- c("EvenJ_site", "EvenJ_Est", "EvenJ_SE", "EvenJ_p")
    colnames(ev_df) <- xn
head(ev_df)

###############################################################
#abundance
hist(log10(abundance))

#subset data for loop
df<-data.frame(site_id, year.1, log10(abundance))
x <- c("site", "yr", "abun")
    colnames(df) <- x
head(df)

df2 <- df[which(df$site!=108000043& df$site!=108000044& 
df$site!=108000057& df$site!=109000270& df$site!=121000206),]

trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(abun ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

df108000043 <- df[which(df$site==108000043),]
g1<-gls(abun ~ yr,na.action=na.omit, data = df108000043)
trend.108000043 <- summary(gls(abun ~ yr,na.action=na.omit, data = df108000043))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df108000043  <- data.frame(site = 108000043, t(trend.108000043))

df108000044 <- df[which(df$site==108000044),]
g1<-gls(abun ~ yr,na.action=na.omit, data = df108000044)
trend.108000044 <- summary(gls(abun ~ yr,na.action=na.omit, data = df108000044))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df108000044 <- data.frame(site = 108000044, t(trend.108000044))

df108000057 <- df[which(df$site==108000057),]
g1<-gls(abun ~ yr,na.action=na.omit, data = df108000057)
trend.108000057 <- summary(gls(abun ~ yr,na.action=na.omit, data = df108000057))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df108000057 <- data.frame(site = 108000057, t(trend.108000057))

df109000270 <- df[which(df$site==109000270),]
g1<-gls(abun ~ yr,na.action=na.omit, data = df109000270)
trend.109000270 <- summary(gls(abun ~ yr,na.action=na.omit, data = df109000270))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df109000270  <- data.frame(site = 109000270, t(trend.109000270))

df121000206 <- df[which(df$site==121000206),]
g1<-gls(abun ~ yr,na.action=na.omit, data = df121000206)
trend.121000206 <- summary(gls(abun ~ yr,na.action=na.omit, data = df121000206))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df121000206  <- data.frame(site = 121000206, t(trend.121000206))

tot_t<-rbind(trends, df108000043,df108000044, df108000057, df109000270, df121000206)
ab_df <- tot_t[order(tot_t$site),] 

xn <- c("Adun_site", "Adun_Est", "Adun_SE", "Adun_p")
    colnames(ab_df) <- xn
nrow(ab_df)
##############################################
#turnover
hist(turnover)

#subset data for loop
df<-data.frame(site_id, year.1, turnover)
x <- c("site", "yr", "to")
    colnames(df) <- x
head(df)

df2 <- df[which(df$site!=121000013),]

trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(to ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

df121000013 <- df[which(df$site==121000013),]
g1<-gls(to ~ yr,na.action=na.omit, data = df121000013)
trend.121000013 <- summary(gls(to ~ yr,na.action=na.omit, data = df121000013))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df121000013  <- data.frame(site = 121000013, t(trend.121000013))

tot_t<-rbind(trends, df121000013)
to_df <- tot_t[order(tot_t$site),] 

xn <- c("TurnO_site", "TurnO_Est", "TurnO_SE", "TurnO_p")
    colnames(to_df) <- xn
nrow(to_df)

##############################################
#F_to
hist(F_to)

#subset data for loop
df<-data.frame(site_id, year.1, F_to)
x <- c("site", "yr", "F_to")
    colnames(df) <- x
head(df)

df2 <- df[which(df$site!=103000552& df$site!=104000141& 
df$site!=104000142& df$site!=108000153),]

trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(F_to ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

df103000552 <- df[which(df$site==103000552),]
g1<-gls(F_to ~ yr,na.action=na.omit, data = df103000552)
#trend.103000552 <- summary(gls(F_to ~ yr,na.action=na.omit, data = df103000552))$tTable[2, c(1,2,4)]
#acf(residuals(g1,type="p")) # no temporal autocorrelation
df103000552 <- data.frame(site = 103000552, Value = 0, Std.Error = "NA", p.value = "NA")

df104000141 <- df[which(df$site==104000141),]
g1<-gls(F_to ~ yr,na.action=na.omit, data = df104000141)
trend.104000141 <- summary(gls(F_to ~ yr,na.action=na.omit, data = df104000141))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df104000141 <- data.frame(site = 104000141, t(trend.104000141))

df104000142 <- df[which(df$site==104000142),]
g1<-gls(F_to ~ yr,na.action=na.omit, data = df104000142)
#trend.104000142 <- summary(gls(F_to ~ yr,na.action=na.omit, data = df104000142))$tTable[2, c(1,2,4)]
#acf(residuals(g1,type="p")) #no temporal autocorrelation
df104000142 <- data.frame(site = 104000142, Value = 0, Std.Error = "NA", p.value = "NA")

df108000153 <- df[which(df$site==108000153),]
g1<-gls(F_to ~ yr,na.action=na.omit, data = df108000153)
trend.108000153 <- summary(gls(F_to ~ yr,na.action=na.omit, data = df108000153))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df108000153  <- data.frame(site = 108000153, t(trend.108000153))

tot_t<-rbind(trends, df103000552,df104000141, df104000142, df108000153)
F_to_df <- tot_t[order(tot_t$site),] 

xn <- c("F_to_site", "F_to_Est", "F_to_SE", "F_to_p")
    colnames(F_to_df) <- xn
nrow(F_to_df)
##############################################
#FRic
hist(FRic)
FR<-(FRic/78543.77) #divide by maximum value to scale

#subset data for loop
df<-data.frame(site_id, year.1, (FRic/78543.77))
df<-data.frame(site_id, year.1, FRic)
x <- c("site", "yr", "FR")
    colnames(df) <- x
head(df)

df2 <- df[which(df$site!=103000103),] # for scaled
##df2 <- df[which(df$site!=100000094),] # for not scaled

trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(FR ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

df103000103 <- df[which(df$site==103000103),]
g1<-gls(FR ~ yr,na.action=na.omit, data = df103000103)
trend.103000103 <- summary(gls(FR ~ yr,na.action=na.omit, data = df103000103))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df103000103  <- data.frame(site = 103000103, t(trend.103000103))

tot_t<-rbind(trends, df103000103)
FRic_df <- tot_t[order(tot_t$site),] 

xn <- c("FRic_site", "FRic_Est", "FRic_SE", "FRic_p")
    colnames(FRic_df) <- xn
nrow(FRic_df)

##############################################
#FEve
hist(FEve)

#subset data for loop
df<-data.frame(site_id, year.1, FEve)
x <- c("site", "yr", "FEve")
    colnames(df) <- x
head(df)

df2 <- df[which(df$site!=103000056& df$site!=114000047),]

trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(FEve ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

df103000056 <- df[which(df$site==103000056),]
g1<-gls(FEve ~ yr,na.action=na.omit, data = df103000056)
trend.103000056 <- summary(gls(FEve ~ yr,na.action=na.omit, data = df103000056))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df103000056 <- data.frame(site = 103000056, t(trend.103000056))

df114000047 <- df[which(df$site==114000047),]
g1<-gls(FEve ~ yr,na.action=na.omit, data = df114000047)
trend.114000047 <- summary(gls(FEve ~ yr,na.action=na.omit, data = df114000047))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df114000047  <- data.frame(site = 114000047, t(trend.114000047))

tot_t<-rbind(trends, df103000056,df114000047)
FEve_df <- tot_t[order(tot_t$site),] 

xn <- c("FEve_site", "FEve_Est", "FEve_SE", "FEve_p")
    colnames(FEve_df) <- xn
nrow(FEve_df)

##############################################
#FDiv
hist(FDiv)

#subset data for loop
df<-data.frame(site_id, year.1, FDiv)
x <- c("site", "yr", "FDiv")
    colnames(df) <- x
head(df)

trends <- NULL
for(i in unique(df$site)){
  sub <- df[df$site == i, ]
  trend.i <- summary(gls(FDiv ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

FDiv_df <- trends[order(trends$site),] 
xn <- c("FDiv_site", "FDiv_Est", "FDiv_SE", "FDiv_p")
    colnames(FDiv_df) <- xn
nrow(FDiv_df)

##############################################
#FDis
hist(FDis)

#subset data for loop
df<-data.frame(site_id, year.1, FDis)
x <- c("site", "yr", "FDis")
    colnames(df) <- x
head(df)

df2 <- df[which(df$site!=100000167& df$site!=108000024),]

trends <- NULL
for(i in unique(df2$site)){
  sub <- df2[df2$site == i, ]
  trend.i <- summary(gls(FDis ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

df100000167 <- df[which(df$site==100000167),]
g1<-gls(FDis ~ yr,na.action=na.omit, data = df100000167)
trend.100000167 <- summary(gls(FDis ~ yr,na.action=na.omit, data = df100000167))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) #no temporal autocorrelation
df100000167 <- data.frame(site = 100000167, t(trend.100000167))

df108000024 <- df[which(df$site==108000024),]
g1<-gls(FDis ~ yr,na.action=na.omit, data = df108000024)
trend.108000024 <- summary(gls(FDis ~ yr,na.action=na.omit, data = df108000024))$tTable[2, c(1,2,4)]
acf(residuals(g1,type="p")) # no temporal autocorrelation
df108000024  <- data.frame(site = 108000024, t(trend.108000024))

tot_t<-rbind(trends, df100000167,df108000024)
FDis_df <- tot_t[order(tot_t$site),] 

xn <- c("FDis_site", "FDis_Est", "FDis_SE", "FDis_p")
    colnames(FDis_df) <- xn
nrow(FDis_df)

##############################################
#RaoQ
hist(RaoQ)

#subset data for loop
df<-data.frame(site_id, year.1, RaoQ)
x <- c("site", "yr", "RaoQ")
    colnames(df) <- x
head(df)

trends <- NULL
for(i in unique(df$site)){
  sub <- df[df$site == i, ]
  trend.i <- summary(gls(RaoQ ~ yr,correlation = corAR1(form = ~ yr),na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
  trend.i <- data.frame(site = i, 
                        t(trend.i))
  trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
} ; rm(i)

RaoQ_df <- trends[order(trends$site),] 
xn <- c("RaoQ_site", "RaoQ_Est", "RaoQ_SE", "RaoQ_p")
    colnames(RaoQ_df) <- xn
nrow(RaoQ_df)
##############################################

glm.output<- cbind(sr_df, si_df,sh_df,ev_df,ab_df,to_df,F_to_df,FRic_df,FEve_df,FDiv_df,RaoQ_df)
head(glm.output)

glm.output$Adun_site<-as.numeric(glm.output$Adun_site)
glm.output$ShanH_site<-as.numeric(glm.output$ShanH_site)
glm.output$EvenJ_site<-as.numeric(glm.output$EvenJ_site)
glm.output$FDiv_site<-as.numeric(glm.output$FDiv_site)
glm.output$RaoQ_site<-as.numeric(glm.output$RaoQ_site)

identical(glm.output$SppRich_site,glm.output$SimpD_site)
identical(glm.output$SppRich_site,glm.output$ShanH_site)
identical(glm.output$SppRich_site,glm.output$EvenJ_site)
identical(glm.output$SppRich_site,glm.output$Adun_site)
identical(glm.output$SppRich_site,glm.output$TurnO_site)
identical(glm.output$SppRich_site,glm.output$F_to_site)
identical(glm.output$SppRich_site,glm.output$FRic_site)
identical(glm.output$SppRich_site,glm.output$FEve_site)
identical(glm.output$SppRich_site,glm.output$FDiv_site)
identical(glm.output$SppRich_site,glm.output$RaoQ_site)

library(dplyr)
new_df<-select(glm.output, -c("SimpD_site","ShanH_site","EvenJ_site","Adun_site",
"TurnO_site","F_to_site","FRic_site","FEve_site","FDiv_site","RaoQ_site"))

names(new_df)[names(new_df) == 'SppRich_site'] <- 'site'
head(new_df)

##
write.csv(new_df,"glmOutput.csv")
##