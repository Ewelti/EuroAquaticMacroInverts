########################################################################################################################################################################
##Set working directory
#setwd("C:/Users/ewelti/Desktop/aquatic_data/data_csvs")

# Load libraries: 
library(vegan)
library(codyn)
library(reshape2)
##########################################
######################################

# Read data:

DATA1_list <- read.csv("raw-data/All_000.csv", header=T) # simple version of all time series abundance site-yr-taxa data
head(DATA1_list)
unique(DATA1_list$site_id)
#
sub_data<-DATA1_list[ which(DATA1_list$site_id=='109000001'),] # change file name according to the analysed time series

# (1) Compute biodiversity metrics ---------------------------------------------------------------
DATA1 <- dcast(sub_data, year ~ taxon_id, sum, value.var = "abundance") # create cross table; yr by spp matrix
lastTaxon <- length(DATA1) #assign the number position of last listed spp to "lastTaxon"
firstTaxon <- 2 #first listed spp is in the second column (after yr)
YR <-DATA1[,1]

DATA1_Taxa <- subset(DATA1[,c(firstTaxon:lastTaxon)]) #creates a sub dataset without yr column 
DATA1$NTaxa <- specnumber(DATA1_Taxa)  # taxonomic richness per yr And adds this col to DATA1
DATA1$Simp <- diversity(DATA1_Taxa, index = "simpson") # Simpson?s taxonomic diversity per yr And adds this col to DATA1
DATA1$Shan <- diversity(DATA1_Taxa, index = "shannon")
DATA1$Even <- DATA1$Shan/log(DATA1$NTaxa) #Pielou's evenness (J)
DATA1$Abund <- rowSums (DATA1_Taxa) # Total abundance per yr And adds this col to DATA1
DATA1_Turnover <- turnover(sub_data, time.var = "year", species.var = "taxon_id", abundance.var = "abundance" , metric = "total")
DATA1$Turnover <- c("NA", DATA1_Turnover$total) # Turnover per yr And adds this col to DATA1

bd_ <- data.frame(year=YR, spp_rich=DATA1$NTaxa,simpD=DATA1$Simp, shanH=DATA1$Shan, evenJ=DATA1$Even, abun=DATA1$Abund, turnover=DATA1$Turnover)
write.csv(bd_, file = "109_001.csv") # change file name according to the analysed time series

#########################################
############################################