###################################################################
#set wd
#setwd("C:/Users/ewelti/Desktop/aquatic_data/data_csvs")
###################################################################

#libraries
library(vegan)
library(reshape2)
library(codyn)
library(mobr)

###############################################################################################################

#load data in longform
all <- read.csv("All_000_noDups.csv", h = T, sep = ",", stringsAsFactors = FALSE, check.names=FALSE) ##this file is on OSF
head(all)

#individual-based rarefication of species richness
all$ro.ab <- round(all$abundance, digits = 0) #rarefication only works on integers

TD <- NULL
for(i in unique(all$site_id)){
  sub <- all[all$site_id == i, ]
  sub.m <- dcast(sub, code ~ taxon_id, sum, value.var = "abundance") 	# matrix form
  sub.ta <- subset(sub.m[,c(2:length(sub.m))])                  		  # subset matrix to remove row names
  SppRich <- specnumber(sub.ta)  							          # taxonomic richness
  Simp <- diversity(sub.ta, index = "simpson") 					# Simpson's taxonomic diversity
  Shan <- diversity(sub.ta, index = "shannon")					# Shannon's taxonomic diversity
  EvenJ <- Shan/log(SppRich) 							    	        # Pielou's evenness (J)
  E10 <- Shan/SppRich 								  	              # Shannon's evenness (E10)
  Abund <- rowSums (sub.ta) 								            # Total abundance
  S_PIE <- calc_PIE(sub.ta, ENS = TRUE)						      # effective number of common species
  DATA1_Turnover <- turnover(sub, time.var = "year", species.var = "taxon_id", abundance.var = "abundance" , metric = "total")
  Turnover <- c("NA", DATA1_Turnover$total) 					  # Turnover per yr And first yr is "NA"
  sub.m_r <- dcast(sub, code ~ taxon_id, sum, value.var = "ro.ab")      # matrix form for rarefaction with rounded richness
  sub.ta_r <- subset(sub.m_r[,c(2:length(sub.m_r))])                  	# subset matrix to remove row names for rounded sppRich    	
  rare.sppRich <- if (min(rowSums(sub.ta_r)) > 10) {
    rarefy(sub.ta_r, sample = min(rowSums(sub.ta_r)))			# rarefy based on min abundance
  } else {rarefy(sub.ta_r, sample = 10)} 					# rarefy based on abund =10 if min is less
  TD.i <- data.frame(sub.m$code, SppRich, Simp, Shan, EvenJ, E10, Abund, S_PIE, Turnover, rare.sppRich)
  TD <- rbind(TD, TD.i) ; rm(TD.i, sub.m, sub.ta, sub, SppRich, Simp, Shan, 
                             EvenJ, E10, Abund, S_PIE, DATA1_Turnover, Turnover, sub.m_r, sub.ta_r, rare.sppRich)
} ; rm(i)

####################################################
write.csv(TD, "TaxaDiversity_siteYr.csv")
###################################################################################

