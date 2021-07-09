#set wd
#setwd("C:/Users/ewelti/Desktop/aquatic_data/")

#############
#load lbries
library(BAT)
library(FD)
library(vegan)
library(codyn)
library(DataCombine)
##### LOAD DATA###########
### COMMUNITY DATA###

# species by trait matrix
comm_traits <- read.csv("outputs/FunctionalTraits/Traits_fuzzy_BiolOnly_85percCoverage.csv", h = T, sep = ",", row.names = 1, stringsAsFactors = FALSE)
colnames(comm_traits)

#Site by species matrix
comm_abund <- read.csv("outputs/FunctionalTraits/mform_abundances.csv", h = T, sep = ",", row.names = 1, stringsAsFactors = FALSE, check.names=FALSE)
head(comm_abund)

#Site Year factors
####NOTE: this file must be in order by year (smallest to largest) within site
facto <- read.csv("outputs/FunctionalTraits/SiteYears_factors.csv", h = T, sep = ",", stringsAsFactors = FALSE)
head(facto)

#traits
cwm_l <- read.csv("outputs/FunctionalTraits/CWMs_longform.csv", h = T, sep = ",", stringsAsFactors = FALSE)
head(cwm_l)
####################################################################################################
####################################################################################################
###############################################################
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="103000662",]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait,] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$code %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long,] #subset by this list

# Create site subset of factor
spp_fact = (facto$code %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact,] #subset by this list

#check if names are the same
identical(row.names(traits_site),colnames(site_s))

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)

#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
New<-c("NA","2003")
s1_tot_FuncTu_ <- InsertRow(s1_tot_FuncTu, NewRow = New, RowNum = 8)
fact_site$F_to <- c("NA", s1_tot_FuncTu_$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
s1_app_FuncTu_ <- InsertRow(s1_app_FuncTu, NewRow = New, RowNum = 8)
fact_site$F_to_app <- c("NA", s1_app_FuncTu_$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
s1_dis_FuncTu_ <- InsertRow(s1_dis_FuncTu, NewRow = New, RowNum = 8)
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu_$disappearance) 
##

tur_o<-cbind(fact_site,FC_biol_ord)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="103000694",]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait,] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$code %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long,] #subset by this list

# Create site subset of factor
spp_fact = (facto$code %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact,] #subset by this list

#check if names are the same
identical(row.names(traits_site),colnames(site_s))

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
New<-c("NA","1996")
s1_tot_FuncTu_ <- InsertRow(s1_tot_FuncTu, NewRow = New, RowNum = 1)
fact_site$F_to <- c("NA", s1_tot_FuncTu_$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
s1_app_FuncTu_ <- InsertRow(s1_app_FuncTu, NewRow = New, RowNum = 1)
fact_site$F_to_app <- c("NA", s1_app_FuncTu_$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
s1_dis_FuncTu_ <- InsertRow(s1_dis_FuncTu, NewRow = New, RowNum = 1)
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu_$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="109000220",]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait,] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$code %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long,] #subset by this list

# Create site subset of factor
spp_fact = (facto$code %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact,] #subset by this list

#check if names are the same
identical(row.names(traits_site),colnames(site_s))

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
New<-c("NA","1996")
s1_tot_FuncTu_ <- InsertRow(s1_tot_FuncTu, NewRow = New, RowNum = 1)
fact_site$F_to <- c("NA", s1_tot_FuncTu_$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
s1_app_FuncTu_ <- InsertRow(s1_app_FuncTu, NewRow = New, RowNum = 1)
fact_site$F_to_app <- c("NA", s1_app_FuncTu_$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
s1_dis_FuncTu_ <- InsertRow(s1_dis_FuncTu, NewRow = New, RowNum = 1)
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu_$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
write.csv(tur_o, "RowTurnoverProbl_sites.csv")
