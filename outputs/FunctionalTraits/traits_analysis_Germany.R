#set wd
#setwd("C:/Users/ewelti/Desktop/aquatic_data/TraitCalculation")

#############
#load lbries
library(BAT)
library(FD)
library(vegan)
library(codyn)

##### LOAD DATA###########
### COMMUNITY DATA###

# species by trait matrix
comm_traits <- read.csv("outputs/FunctionalTraits/Traits_fuzzy_BiolOnly_85percCoverage.csv", h = T, sep = ",", row.names = 1, stringsAsFactors = FALSE)
colnames(comm_traits)

####NOTE: all these files must be in order by year (smallest to largest) within site##############################
#Site by species matrix
comm_abund <- read.csv("outputs/FunctionalTraits/mform_abundances.csv", h = T, sep = ",", row.names = 1, stringsAsFactors = FALSE, check.names=FALSE)
head(comm_abund)

#Site Year factors
facto <- read.csv("outputs/FunctionalTraits/SiteYears_factors.csv", h = T, sep = ",", stringsAsFactors = FALSE)
head(facto)

#traits long form
cwm_l <- read.csv("outputs/FunctionalTraits/CWMs_longform.csv", h = T, sep = ",", stringsAsFactors = FALSE)
head(cwm_l)
####################################################################################################
####################################################################################################
###############################################################
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000001",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##
tur_o<-cbind(fact_site,FC_biol_ord)

###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000002",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000003",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000004",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
##
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000005",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)

###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000006",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000007",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
##
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000008",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
##
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000009",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)

###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000010",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
##
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000011",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)

###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000012",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000013",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000014",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
##
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000015",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000016",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000017",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000019",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000020",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000021",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000022",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000024",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000025",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000026",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000027",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000028",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000029",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000030",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000031",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000032",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000033",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000034",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000035",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000036",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000037",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000038",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000039",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000040",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000041",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000042",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000043",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000044",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000045",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000046",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000047",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000049",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000050",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000051",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000052",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000053",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000054",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000056",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000057",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000058",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000059",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000060",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000061",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000062",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000063",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000064",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000065",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000066",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000067",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000068",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000069",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000070",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000071",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000072",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000073",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000074",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000075",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000076",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000077",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000078",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000079",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000080",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000081",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000082",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000083",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000084",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000085",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000086",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000087",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000088",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000089",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000090",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000091",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000092",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000093",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000094",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000095",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000096",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000097",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000098",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000099",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000100",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000101",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000102",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000103",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000104",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000106",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000107",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
##
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000108",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
##
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000109",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)

###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000110",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
##
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000111",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)

###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000112",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000113",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000114",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
##
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000115",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000116",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000117",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000118",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000119",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000120",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000121",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000122",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000123",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000124",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000125",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000126",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000127",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000128",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000129",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000130",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000131",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000132",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000133",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000134",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000135",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000136",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000137",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000138",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000139",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000142",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000143",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000144",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000145",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000146",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000147",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000148",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000149",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000150",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000151",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000152",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000153",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000156",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000159",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000162",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000164",]
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

#FunctionalTraits
#principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site,site_s, calc.FRic =T, m = 6, stand.FRic = F, 
                    scale.RaoQ = F, calc.FDiv = T,asym.bin = 5:10,calc.CWM=FALSE)
#
#Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long,time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###########################################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="108000165",]
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
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality",abundance.var = "CWM",metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 
##

tur<-cbind(fact_site,FC_biol_ord)
tur_o<-rbind(tur,tur_o)
###############################################################
write.csv(tur_o, "Germany_FD.csv")