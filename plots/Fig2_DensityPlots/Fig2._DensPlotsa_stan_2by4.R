##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

# attach data
response_stan <- readRDS("outputs/outputs_metaAnalysis/stanTrends_site_level.rds")
head(response_stan)
unique(response_stan$Response)
length(unique(response_stan$site_id))

metricdata <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv")
head(metricdata)

Ests <- read.csv("outputs/outputs_metaAnalysis/Yr_metaanaly_Ests.csv")

#tiff(filename = "plots/Fig2_DensityPlots/Fig2._SlopeDistributions.tiff", width = 10, height = 6, units = 'in', res = 600, compression = 'lzw')
pdf("plots/Fig2_DensityPlots/Fig2._SlopeDistributions.pdf", width = 10, height = 6)

par(mar=c(4,0.4,0.4,0.4), mfrow=c(1,2))

#####################Taxonomic metrics ##############################
#### Spp Richness #####
SppRich <- subset(response_stan, Response == "spp_richness")
SppRich <- SppRich$estimate[!is.na(SppRich$estimate)]
ave_SppRich <- mean(metricdata$spp_richness, na.rm=T)
percChange_perYr<-(SppRich/ave_SppRich)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
#axis(2, at=0, labels="SR", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
SppRich_Est <- subset(Ests, Response == "spp_richness")
stand_SppRich <- lapply(SppRich_Est[,2:11],"*",100/ave_SppRich)
yy <- (3.2/4*(b-a)+a)
points(x=stand_SppRich$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_SppRich$Q2.5, stand_SppRich$Q2.5, stand_SppRich$Q97.5, stand_SppRich$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_SppRich$Q5, stand_SppRich$Q5, stand_SppRich$Q95, stand_SppRich$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_SppRich$Q10, stand_SppRich$Q10, stand_SppRich$Q90, stand_SppRich$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd = 1)
legend(x=-24, y=(4/4*(b-a)+a), legend=("a, Taxon richness"), bty="n", cex=1.3)
legend(x=5, y=(4/4*(b-a)+a), legend=expression(paste("+0.73% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Abundance #####
abund <- subset(response_stan, Response == "abundance")
abund <- abund$estimate[!is.na(abund$estimate)]
percChange_perYr <- (10^abund-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
abund_Est <- subset(Ests, Response == "abundance")
stand_abund <- lapply(10^abund_Est[,2:11]-1,"*",100)
yy <- (2.2/4*(b-a)+a)
points(x=stand_abund$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_abund$Q2.5, stand_abund$Q2.5, stand_abund$Q97.5, stand_abund$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_abund$Q5, stand_abund$Q5, stand_abund$Q95, stand_abund$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_abund$Q10, stand_abund$Q10, stand_abund$Q90, stand_abund$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(3/4*(b-a)+a), legend=("b, Abundance"), bty="n", cex=1.3)
legend(x=5, y=(3/4*(b-a)+a), legend=expression(paste("+1.17% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Shannon's Evenness #####
E10 <- subset(response_stan, Response == "E10")
E10 <- E10$estimate[!is.na(E10$estimate)]
percChange_perYr<-(10^E10-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*-1
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="Even", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
E10_Est <- subset(Ests, Response == "E10")
stand_E10 <- lapply(10^E10_Est[,2:11]-1,"*",100)
yy <- (1.2/4*(b-a)+a)
points(x=stand_E10$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_E10$Q2.5, stand_E10$Q2.5, stand_E10$Q97.5, stand_E10$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_E10$Q5, stand_E10$Q5, stand_E10$Q95, stand_E10$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_E10$Q10, stand_E10$Q10, stand_E10$Q90, stand_E10$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(2/4*(b-a)+a), legend=("c, Evenness"), bty="n", cex=1.3)
legend(x=5, y=(2/4*(b-a)+a), legend=expression(paste("-0.54% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Turnover #####
turnover <- subset(response_stan, Response == "turnover")
turnover <- turnover$estimate[!is.na(turnover$estimate)]
to<-as.numeric(replace(metricdata$turnover, metricdata$turnover=="<NA>", "NA"))
ave_turnover <- mean(to, na.rm=T)
percChange_perYr<-(turnover/ave_turnover)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="TurnO", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
turnover_Est <- subset(Ests, Response == "turnover")
stand_turnover <- lapply(turnover_Est[,2:11],"*",100/ave_turnover)
yy <- (0.2/4*(b-a)+a)
points(x=stand_turnover$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_turnover$Q2.5, stand_turnover$Q2.5, stand_turnover$Q97.5, stand_turnover$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_turnover$Q5, stand_turnover$Q5, stand_turnover$Q95, stand_turnover$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_turnover$Q10, stand_turnover$Q10, stand_turnover$Q90, stand_turnover$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(1/4*(b-a)+a), legend=("d, Turnover"), bty="n", cex=1.3)
legend(x=5, y=(1/4*(b-a)+a), legend=expression(paste("-0.2% ",y^-1,, sep = "")), bty="n", cex=1.3)

##
box(lwd=2,col="white")
axis(1,lwd=2)
abline(v=0, lwd=1.5, lty=2)
##

#####################Functional metrics ##############################
#### Functional Richness #####
FRic <- subset(response_stan, Response == "FRic")
FRic <- FRic$estimate[!is.na(FRic$estimate)]
percChange_perYr<-(10^FRic-1)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',yaxs="i")
title(xlab=expression(paste("% change y"^"-1")), line=2.4,cex.lab=1.3)
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
FRic_Est <- subset(Ests, Response == "func_rich")
stand_FRic <- lapply(10^FRic_Est[,2:11]-1,"*", 100)
yy <- (3.2/4*(b-a)+a)
points(x=stand_FRic$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FRic$Q2.5, stand_FRic$Q2.5, stand_FRic$Q97.5, stand_FRic$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FRic$Q5, stand_FRic$Q5, stand_FRic$Q95, stand_FRic$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FRic$Q10, stand_FRic$Q10, stand_FRic$Q90, stand_FRic$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(4/4*(b-a)+a), legend=("e, Func. richness"), bty="n", cex=1.3)
legend(x=5, y=(4/4*(b-a)+a), legend=expression(paste("+2.4% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Functional Redundancy #####
FRed <- subset(response_stan, Response == "FRed")
FRed <- FRed$estimate[!is.na(FRed$estimate)]
ave_FRed <- mean(metricdata$FRed, na.rm=T)
percChange_perYr<-(FRed/ave_FRed)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
FRed_Est <- subset(Ests, Response == "FRed")
stand_FRed <- lapply(FRed_Est[,2:11],"*",100/ave_FRed)
yy <- (2.2/4*(b-a)+a)
points(x=stand_FRed$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FRed$Q2.5, stand_FRed$Q2.5, stand_FRed$Q97.5, stand_FRed$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FRed$Q5, stand_FRed$Q5, stand_FRed$Q95, stand_FRed$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FRed$Q10, stand_FRed$Q10, stand_FRed$Q90, stand_FRed$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(3/4*(b-a)+a), legend=("f, Func. redundancy"), bty="n", cex=1.3)
legend(x=5, y=(3/4*(b-a)+a), legend=expression(paste("+0.03% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Func evenness #####
FEve <- subset(response_stan, Response == "FEve")
FEve <- FEve$estimate[!is.na(FEve$estimate)]
ave_FEve <- mean(metricdata$FEve, na.rm=T)
percChange_perYr<-(FEve/ave_FEve)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-1
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
#axis(2, at=0, labels="FEve", las=1,cex.axis=1.3)
##
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =2)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)
FEve_Est <- subset(Ests, Response == "func_even")
stand_FEve <- lapply(FEve_Est[,2:11],"*",100/ave_FEve)
yy <- (1.2/4*(b-a)+a)
points(x=stand_FEve$Estimate, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stand_FEve$Q2.5, stand_FEve$Q2.5, stand_FEve$Q97.5, stand_FEve$Q97.5),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stand_FEve$Q5, stand_FEve$Q5, stand_FEve$Q95, stand_FEve$Q95),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stand_FEve$Q10, stand_FEve$Q10, stand_FEve$Q90, stand_FEve$Q90),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(2/4*(b-a)+a), legend=("g, Func. evenness"), bty="n", cex=1.3)
legend(x=5, y=(2/4*(b-a)+a), legend=expression(paste("-0.22% ",y^-1,, sep = "")), bty="n", cex=1.3)

#### Func turnover #####
response_fto <- readRDS("outputs/meta_updated/metaanalysis_F_to.rds")
F_to <- subset(response_stan, Response == "F_to")

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
ft_merge <- merge(F_to,ft_merg,by="site_id")
head(ft_merge)

##predict % change/yr for each sites
library(betareg)

percCh <- NULL
for(i in unique(ft_merge$site_id)){
    sub <- ft_merge[ft_merge$site_id == i, ]
	yr.1 <- sub$med_yr- sub$max_yr
	yr.n <- sub$med_yr - sub$min_yr
	year1 <- plogis(sub$intercept + sub$estimate * yr.1)
	yearn <- plogis(sub$intercept + sub$estimate * yr.n)
	rate <- (yearn/year1)^(1/(test$max_yr-test$min_yr))
	perc = (rate -1)*100
    percCh.i <- data.frame(site_id = i, perc)
    percCh <- rbind(percCh, percCh.i) ; rm(percCh.i, sub, year1, yearn, rate, perc)
} ; rm(i)
head(percCh)

percChange_perYr <- na.omit(percCh$perc)
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *0
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-20,20),ylim=c(a,b),col="white",yaxt='n',xaxt='n',yaxs="i")
polygon(c(d$x[d$x >= 0 ], 0),
        c(d$y[d$x >= 0 ], 0),
        col = "cornflowerblue", border = "cornflowerblue", lwd =3)
polygon(c(d$x[d$x <= 0 ], 0),
        c(d$y[d$x <= 0 ], 0),
        col = "coral1", border = "coral1", lwd =2)


F_to_Est <- subset(Ests, Response == "func_turnover")
head(F_to_Est)

allinter <- mean(ft_merge$intercept)
minyr <- min(ft_merge$min_yr)
maxyr <- max(ft_merge$max_yr)

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

stand_F_to <- stan_Fto$perc[stan_Fto$output=="Estimate"]
yy <- (0.2/4*(b-a)+a)
points(x=stand_F_to, y=yy, lwd=2,pch="|",cex=2)
polygon(x=c(stan_Fto$perc[stan_Fto$output=="Q2.5"], stan_Fto$perc[stan_Fto$output=="Q2.5"], stan_Fto$perc[stan_Fto$output=="Q97.5"], stan_Fto$perc[stan_Fto$output=="Q97.5"]),
        y=c((yy-yy/18),(yy+yy/18),(yy+yy/18),(yy-yy/18)),
        col = 1,border = 0,lwd = 1)
polygon(x=c(stan_Fto$perc[stan_Fto$output=="Q5"], stan_Fto$perc[stan_Fto$output=="Q5"], stan_Fto$perc[stan_Fto$output=="Q95"], stan_Fto$perc[stan_Fto$output=="Q95"]),
        y=c((yy-yy/10),(yy+yy/10),(yy+yy/10),(yy-yy/10)),
        col = 1,border = 0,lwd =1)
polygon(x=c(stan_Fto$perc[stan_Fto$output=="Q10"], stan_Fto$perc[stan_Fto$output=="Q10"], stan_Fto$perc[stan_Fto$output=="Q90"], stan_Fto$perc[stan_Fto$output=="Q90"]),
        y=c((yy-yy/6),(yy+yy/6),(yy+yy/6),(yy-yy/6)),
        col = 1,border = 0,lwd =1)
legend(x=-24, y=(1/4*(b-a)+a), legend=("h, Func. turnover"), bty="n", cex=1.3)
legend(x=5, y=(1/4*(b-a)+a), legend=expression(paste("-0.32% ",y^-1,, sep = "")), bty="n", cex=1.3)

##
box(lwd=2,col="white")
axis(1,lwd=2)
abline(v=0, lwd=1.5, lty=2)
##

dev.off()
########################################