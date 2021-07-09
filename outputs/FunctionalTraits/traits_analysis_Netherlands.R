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
site_sub <- comm_abund[facto$site_id=="117000001",]
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
##############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="117000002",]
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
site_sub <- comm_abund[facto$site_id=="117000003",]
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
site_sub <- comm_abund[facto$site_id=="117000004",]
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
site_sub <- comm_abund[facto$site_id=="117000005",]
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
site_sub <- comm_abund[facto$site_id=="117000006",]
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
site_sub <- comm_abund[facto$site_id=="117000007",]
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
site_sub <- comm_abund[facto$site_id=="117000008",]
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
site_sub <- comm_abund[facto$site_id=="117000009",]
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
site_sub <- comm_abund[facto$site_id=="117000010",]
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
site_sub <- comm_abund[facto$site_id=="117000011",]
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
site_sub <- comm_abund[facto$site_id=="117000012",]
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
site_sub <- comm_abund[facto$site_id=="117000013",]
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
site_sub <- comm_abund[facto$site_id=="117000014",]
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
site_sub <- comm_abund[facto$site_id=="117000015",]
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
site_sub <- comm_abund[facto$site_id=="117000016",]
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
site_sub <- comm_abund[facto$site_id=="117000017",]
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
site_sub <- comm_abund[facto$site_id=="117000018",]
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
site_sub <- comm_abund[facto$site_id=="117000019",]
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
site_sub <- comm_abund[facto$site_id=="117000020",]
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
###########################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="117000021",]
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
site_sub <- comm_abund[facto$site_id=="117000022",]
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
site_sub <- comm_abund[facto$site_id=="117000023",]
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
site_sub <- comm_abund[facto$site_id=="117000024",]
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
site_sub <- comm_abund[facto$site_id=="117000025",]
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
site_sub <- comm_abund[facto$site_id=="117000027",]
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
site_sub <- comm_abund[facto$site_id=="117000028",]
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
site_sub <- comm_abund[facto$site_id=="117000029",]
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
site_sub <- comm_abund[facto$site_id=="117000030",]
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
site_sub <- comm_abund[facto$site_id=="117000031",]
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
site_sub <- comm_abund[facto$site_id=="117000032",]
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
site_sub <- comm_abund[facto$site_id=="117000033",]
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
site_sub <- comm_abund[facto$site_id=="117000034",]
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
site_sub <- comm_abund[facto$site_id=="117000035",]
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
site_sub <- comm_abund[facto$site_id=="117000036",]
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
site_sub <- comm_abund[facto$site_id=="117000037",]
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
site_sub <- comm_abund[facto$site_id=="117000038",]
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
site_sub <- comm_abund[facto$site_id=="117000039",]
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
site_sub <- comm_abund[facto$site_id=="117000040",]
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
site_sub <- comm_abund[facto$site_id=="117000041",]
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
site_sub <- comm_abund[facto$site_id=="117000042",]
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
site_sub <- comm_abund[facto$site_id=="117000043",]
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
site_sub <- comm_abund[facto$site_id=="117000044",]
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
site_sub <- comm_abund[facto$site_id=="117000045",]
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
site_sub <- comm_abund[facto$site_id=="117000047",]
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
site_sub <- comm_abund[facto$site_id=="117000048",]
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
site_sub <- comm_abund[facto$site_id=="117000049",]
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
site_sub <- comm_abund[facto$site_id=="117000051",]
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
site_sub <- comm_abund[facto$site_id=="117000052",]
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
site_sub <- comm_abund[facto$site_id=="117000053",]
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
site_sub <- comm_abund[facto$site_id=="117000054",]
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
write.csv(tur_o,"Netherlands_FD.csv")