##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs")

# load libraries
#install.packages("psych")
library(psych)

# attach data
DATA1 <- read.csv("All_siteLevel_and_glmOutput.csv", header=T)
head(DATA1)
DATA1$tmax_est <-as.numeric(DATA1$tmax_est)

cor(DATA1$tmax_Est,DATA1$elevation_atPoint,  method = ("pearson"))

DATA1_list <- read.csv("All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T)
DATA2 <- DATA1_list[!is.na(DATA1_list$site_id_wMissing),]
DATA2$turnover <- as.numeric(DATA2$turnover)
attach(DATA2)
head(DATA2)

corr(

td <- cbind.data.frame(spp_richness,spp_rich_rare,shannonsH, E10,log10(abundance),log10(S_PIE),turnover)
head(td)

pairs.panels(td, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             )

fd <- cbind.data.frame(log10(FRic),FRed,FDiv,RaoQ,FEve,F_to)
head(fd)

pairs.panels(fd, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             )

##all
SR <- spp_richness
RSR <- spp_rich_rare
ShH <- shannonsH
Abun <- log10(abundance)
E10 <- log10(DATA2$E10+0.01)
TO <- turnover
DATA2$FRic <- log10(DATA2$FRic+0.1)
FTO <- F_to
dd <- cbind.data.frame(SR,RSR,ShH, E10, Abun, TO, FRic, FRed, FDiv, RaoQ, FEve, FTO)
head(dd)

dd <- cbind.data.frame(SR, FRic, FRed)
head(dd)

pairs.panels(dd, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             )

##TD and alien/native subsets

alienSR <- alien_SppRich
alienAb <- log10(alien_Abund +1)
natSR <- SppRich_nativeSpp
natAb <- log10(abund_nativeSpp +1)


dd <- cbind.data.frame(SR,RSR,ShH, E10, Abun, TO, alienSR, alienAb, natSR, natAb)
head(dd)

pairs.panels(dd, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
)

##TD and EPT/insect subsets

EPTSR <- EPT_SppRich
EPTAb <- log10(EPT_Abund +1)
insSR <- insect_SppRich
insAb <- log10(insect_Abund +1)


dd <- cbind.data.frame(SR,RSR,ShH, E10, Abun, TO, EPTSR, EPTAb, insSR, insAb)
head(dd)

pairs.panels(dd, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
)

##FD and alien/native subsets

dd <- cbind.data.frame(FRic, FRed, FDiv, RaoQ, FEve, FTO, alienSR, alienAb, natSR, natAb)
head(dd)

pairs.panels(dd, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
)

##FD and EPT/insect subsets

dd <- cbind.data.frame(FRic, FRed, FDiv, RaoQ, FEve, FTO, EPTSR, EPTAb, insSR, insAb)
head(dd)

pairs.panels(dd, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
)

########################################
##response trend correlations
DA <- read.csv("Stan_trends.csv", header=T)
DATA2$turnover <- as.numeric(DATA2$turnover)
attach(DA)
head(DA)

##all
SR <- spp_richness
RSR <- spp_rich_rare
ShH <- shannonsH
Abun <- abundance
TO <- turnover
FTO <- F_to
alienSR <- alien_SppRich
alienAb <- alien_Abund
natSR <- SppRich_nativeSpp
natAb <- abund_nativeSpp
EPTSR <- EPT_SppRich
EPTAb <- EPT_Abund
insSR <- insect_SppRich
insAb <- insect_Abund
dd <- cbind.data.frame(SR,RSR,ShH, E10, Abun, TO, FRic, FRed, FDiv, RaoQ, FEve, FTO, alienSR, alienAb, natSR, natAb, EPTSR, EPTAb, insSR, insAb)
head(dd)

pairs.panels(dd, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = F,  # show density plots
		 ellipses=F,
		 cex.cor =1.5
)

#######################################
##import metadata with sites locations
# attach data
sites <- read.csv("All_siteLevel_and_glmOutput.csv", header=T) # change file name according to the time series to be analyzed
head(sites)
attach(sites)

ff <- lm(tmax_Est~elevation_atPoint)
summary (ff)
plot(tmax_Est~elevation_atPoint)

climate <- cbind.data.frame(aet_Est,ppt_Est,tmax_Est,tmin_Est,aet_mm_12moPrior, ppt_mm_12moPrior, tmax_C_12moPrior, tmin_C_12moPrior)
head(climate)

plot(aet_Est~aet_mm_12moPrior)
pairs.panels(climate, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             )

##
streamOrder <- as.numeric(strahler_streamOrder)
str <- cbind.data.frame(streamOrder,log10(accumulation_atPoint),log10(elevation_atPoint),log10(slope_atPoint+0.001))
head(str)

pairs.panels(str, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             )

##
nitrogen <- cbind.data.frame(NH4_Est,NO3_Est,N_Est,NH4_mean,NO3_mean, N_mean)
head(nitrogen)

pairs.panels(nitrogen, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             )
##
land <- cbind.data.frame(crop_Est, forest_Est, grass_Est, urban_Est,crop_meanPerc_upstr,forest_meanPerc_upstr,grassShrubland_meanPerc_upstr,urban_meanPerc_upstr)
head(land)

pairs.panels(land, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             )
##
dam <- cbind.data.frame(dam_impact_score_lessthan100km,dam_num_connected_lessthan100km, dam_minDist_km_lessthan100km)
head(dam)

pairs.panels(dam, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             )
##
ppt_slope <-ppt_Est
tmax_slope <- tmax_Est
ppt_mean <- ppt_mm_12moPrior
tmax_mean <- tmax_C_12moPrior
strOrder <- strahler_streamOrder
accum <- accumulation_atPoint
elev <- elevation_atPoint
slope <- slope_mean
N_slope <- N_Est
urban_mean <- urban_meanPerc_upstr
crop_mean <- crop_meanPerc_upstr
crop_slope <- crop_Est
urban_slope <- urban_Est
dam_impact <- dam_impact_score_lessthan100km

drivers <- cbind.data.frame(ppt_slope, tmax_slope, ppt_mean, tmax_mean, 
              strOrder, accum, elev, slope, urban_mean, crop_mean, dam_impact) #N_slope, N_mean, crop_slope, urban_slope, 
             


pairs.panels(drivers, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             )
#######################################################################################

all <- cbind.data.frame(SppRich_Est, SppRichRare_Est, ShanH_Est, Abun_Est, TurnO_Est, F_to_Est, FRic_Est, FEve_Est, FDiv_Est, RaoQ_Est,
     		EPT_SppRich_Est, EPT_Abund_Est, insect_SppRich_Est, insect_Abund_Est, AlienSppRich_Est, AlienAbun_Est, nativeSppRich_Est, 
		nativeAbun_Est, ppt_Est, tmax_Est, ppt_mm_12moPrior, tmax_C_12moPrior, 
              strahler_streamOrder, accumulation_atPoint, elevation_atPoint,
              slope_mean, N_Est, N_mean, urban_meanPerc_upstr, crop_meanPerc_upstr,
              crop_Est, urban_Est, dam_impact_score_lessthan100km)

par(mfrow=c(3,5))
# spp rich
plot(all[,1]~all[,19], xlab=colnames(all[19]), ylab="Precip slope")
plot(all[,1]~all[,20], xlab=colnames(all[20]), ylab="Tmax slope")
plot(all[,1]~all[,21], xlab=colnames(all[21]), ylab="Precip mean (mm)")
plot(all[,1]~all[,22], xlab=colnames(all[22]), ylab="Tmax mean (C)")
plot(all[,1]~all[,23], xlab=colnames(all[23]), ylab=")
plot(all[,1]~all[,24], xlab=colnames(all[24]), ylab=colnames(all[1]))
plot(all[,1]~all[,25], xlab=colnames(all[25]), ylab=colnames(all[1]))
plot(all[,1]~all[,26], xlab=colnames(all[26]), ylab=colnames(all[1]))
plot(all[,1]~all[,27], xlab=colnames(all[27]), ylab=colnames(all[1]))
plot(all[,1]~all[,28], xlab=colnames(all[28]), ylab=colnames(all[1]))
plot(all[,1]~all[,29], xlab=colnames(all[29]), ylab=colnames(all[1]))
plot(all[,1]~all[,30], xlab=colnames(all[30]), ylab=colnames(all[1]))
plot(all[,1]~all[,31], xlab=colnames(all[31]), ylab=colnames(all[1]))
plot(all[,1]~all[,32], xlab=colnames(all[32]), ylab=colnames(all[1]))
plot(all[,1]~all[,33], xlab=colnames(all[33]), ylab=colnames(all[1]))

# spp rich rarefied
plot(all[,2]~all[,19], xlab=colnames(all[19]), ylab=colnames(all[2]))
plot(all[,2]~all[,20], xlab=colnames(all[20]), ylab=colnames(all[2]))
plot(all[,2]~all[,21], xlab=colnames(all[21]), ylab=colnames(all[2]))
plot(all[,2]~all[,22], xlab=colnames(all[22]), ylab=colnames(all[2]))
plot(all[,2]~all[,23], xlab=colnames(all[23]), ylab=colnames(all[2]))
plot(all[,2]~all[,24], xlab=colnames(all[24]), ylab=colnames(all[2]))
plot(all[,2]~all[,25], xlab=colnames(all[25]), ylab=colnames(all[2]))
plot(all[,2]~all[,26], xlab=colnames(all[26]), ylab=colnames(all[2]))
plot(all[,2]~all[,27], xlab=colnames(all[27]), ylab=colnames(all[2]))
plot(all[,2]~all[,28], xlab=colnames(all[28]), ylab=colnames(all[2]))
plot(all[,2]~all[,29], xlab=colnames(all[29]), ylab=colnames(all[2]))
plot(all[,2]~all[,30], xlab=colnames(all[30]), ylab=colnames(all[2]))
plot(all[,2]~all[,31], xlab=colnames(all[31]), ylab=colnames(all[2]))
plot(all[,2]~all[,32], xlab=colnames(all[32]), ylab=colnames(all[2]))
plot(all[,2]~all[,33], xlab=colnames(all[33]), ylab=colnames(all[2]))

#Shannon's H
plot(all[,3]~all[,19], xlab=colnames(all[19]), ylab=colnames(all[3]))
plot(all[,3]~all[,20], xlab=colnames(all[20]), ylab=colnames(all[3]))
plot(all[,3]~all[,21], xlab=colnames(all[21]), ylab=colnames(all[3]))
plot(all[,3]~all[,22], xlab=colnames(all[22]), ylab=colnames(all[3]))
plot(all[,3]~all[,23], xlab=colnames(all[23]), ylab=colnames(all[3]))
plot(all[,3]~all[,24], xlab=colnames(all[24]), ylab=colnames(all[3]))
plot(all[,3]~all[,25], xlab=colnames(all[25]), ylab=colnames(all[3]))
plot(all[,3]~all[,26], xlab=colnames(all[26]), ylab=colnames(all[3]))
plot(all[,3]~all[,27], xlab=colnames(all[27]), ylab=colnames(all[3]))
plot(all[,3]~all[,28], xlab=colnames(all[28]), ylab=colnames(all[3]))
plot(all[,3]~all[,29], xlab=colnames(all[29]), ylab=colnames(all[3]))
plot(all[,3]~all[,30], xlab=colnames(all[30]), ylab=colnames(all[3]))
plot(all[,3]~all[,31], xlab=colnames(all[31]), ylab=colnames(all[3]))
plot(all[,3]~all[,32], xlab=colnames(all[32]), ylab=colnames(all[3]))
plot(all[,3]~all[,33], xlab=colnames(all[33]), ylab=colnames(all[3]))

#Abundance
plot(all[,4]~all[,19], xlab=colnames(all[19]), ylab=colnames(all[4]))
plot(all[,4]~all[,20], xlab=colnames(all[20]), ylab=colnames(all[4]))
plot(all[,4]~all[,21], xlab=colnames(all[21]), ylab=colnames(all[4]))
plot(all[,4]~all[,22], xlab=colnames(all[22]), ylab=colnames(all[4]))
plot(all[,4]~all[,23], xlab=colnames(all[23]), ylab=colnames(all[4]))
plot(all[,4]~all[,24], xlab=colnames(all[24]), ylab=colnames(all[4]))
plot(all[,4]~all[,25], xlab=colnames(all[25]), ylab=colnames(all[4]))
plot(all[,4]~all[,26], xlab=colnames(all[26]), ylab=colnames(all[4]))
plot(all[,4]~all[,27], xlab=colnames(all[27]), ylab=colnames(all[4]))
plot(all[,4]~all[,28], xlab=colnames(all[28]), ylab=colnames(all[4]))
plot(all[,4]~all[,29], xlab=colnames(all[29]), ylab=colnames(all[4]))
plot(all[,4]~all[,30], xlab=colnames(all[30]), ylab=colnames(all[4]))
plot(all[,4]~all[,31], xlab=colnames(all[31]), ylab=colnames(all[4]))
plot(all[,4]~all[,32], xlab=colnames(all[32]), ylab=colnames(all[4]))
plot(all[,4]~all[,33], xlab=colnames(all[33]), ylab=colnames(all[4]))

#Turnover
plot(all[,5]~all[,19], xlab=colnames(all[19]), ylab=colnames(all[5]))
plot(all[,5]~all[,20], xlab=colnames(all[20]), ylab=colnames(all[5]))
plot(all[,5]~all[,21], xlab=colnames(all[21]), ylab=colnames(all[5]))
plot(all[,5]~all[,22], xlab=colnames(all[22]), ylab=colnames(all[5]))
plot(all[,5]~all[,23], xlab=colnames(all[23]), ylab=colnames(all[5]))
plot(all[,5]~all[,24], xlab=colnames(all[24]), ylab=colnames(all[5]))
plot(all[,5]~all[,25], xlab=colnames(all[25]), ylab=colnames(all[5]))
plot(all[,5]~all[,26], xlab=colnames(all[26]), ylab=colnames(all[5]))
plot(all[,5]~all[,27], xlab=colnames(all[27]), ylab=colnames(all[5]))
plot(all[,5]~all[,28], xlab=colnames(all[28]), ylab=colnames(all[5]))
plot(all[,5]~all[,29], xlab=colnames(all[29]), ylab=colnames(all[5]))
plot(all[,5]~all[,30], xlab=colnames(all[30]), ylab=colnames(all[5]))
plot(all[,5]~all[,31], xlab=colnames(all[31]), ylab=colnames(all[5]))
plot(all[,5]~all[,32], xlab=colnames(all[32]), ylab=colnames(all[5]))
plot(all[,5]~all[,33], xlab=colnames(all[33]), ylab=colnames(all[5]))

#Functional Richness
plot(all[,6]~all[,19], xlab=colnames(all[19]), ylab=colnames(all[6]))
plot(all[,6]~all[,20], xlab=colnames(all[20]), ylab=colnames(all[6]))
plot(all[,6]~all[,21], xlab=colnames(all[21]), ylab=colnames(all[6]))
plot(all[,6]~all[,22], xlab=colnames(all[22]), ylab=colnames(all[6]))
plot(all[,6]~all[,23], xlab=colnames(all[23]), ylab=colnames(all[6]))
plot(all[,6]~all[,24], xlab=colnames(all[24]), ylab=colnames(all[6]))
plot(all[,6]~all[,25], xlab=colnames(all[25]), ylab=colnames(all[6]))
plot(all[,6]~all[,26], xlab=colnames(all[26]), ylab=colnames(all[6]))
plot(all[,6]~all[,27], xlab=colnames(all[27]), ylab=colnames(all[6]))
plot(all[,6]~all[,28], xlab=colnames(all[28]), ylab=colnames(all[6]))
plot(all[,6]~all[,29], xlab=colnames(all[29]), ylab=colnames(all[6]))
plot(all[,6]~all[,30], xlab=colnames(all[30]), ylab=colnames(all[6]))
plot(all[,6]~all[,31], xlab=colnames(all[31]), ylab=colnames(all[6]))
plot(all[,6]~all[,32], xlab=colnames(all[32]), ylab=colnames(all[6]))
plot(all[,6]~all[,33], xlab=colnames(all[33]), ylab=colnames(all[6]))
