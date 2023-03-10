#set wd
setwd("C:/Users/ewelti/Desktop/aquatic_data/TraitCalculation")

#############
#load lbries
library(BAT)
library(FD)
library(vegan)
library(codyn)

##### LOAD DATA###########
### COMMUNITY DATA###

# species by trait matrix
comm_traits <- read.csv("Traits_fuzzy_BiolOnly_85percCoverage.csv", h = T, sep = ",", row.names = 1, stringsAsFactors = FALSE)
colnames(comm_traits)

####NOTE: all these files must be in order by year (smallest to largest) within site##############################
#Site by species matrix
comm_abund <- read.csv("mform_abundances.csv", h = T, sep = ",", row.names = 1, stringsAsFactors = FALSE, check.names=FALSE)
head(comm_abund)

#Site Year factors
facto <- read.csv("SiteYears_factors.csv", h = T, sep = ",", stringsAsFactors = FALSE)
head(facto)

#traits long form
cwm_l <- read.csv("CWMs_longform.csv", h = T, sep = ",", stringsAsFactors = FALSE)
head(cwm_l)
####################################################################################################
####################################################################################################
###############################################################
###############################################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id=="102000002",]
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
site_sub <- comm_abund[facto$site_id=="102000003",]
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
site_sub <- comm_abund[facto$site_id=="105000001",]
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
site_sub <- comm_abund[facto$site_id=="105000002",]
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
site_sub <- comm_abund[facto$site_id=="105000003",]
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
site_sub <- comm_abund[facto$site_id=="105000004",]
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
site_sub <- comm_abund[facto$site_id=="105000005",]
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
site_sub <- comm_abund[facto$site_id=="105000006",]
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
site_sub <- comm_abund[facto$site_id=="105000007",]
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
site_sub <- comm_abund[facto$site_id=="105000008",]
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
site_sub <- comm_abund[facto$site_id=="105000009",]
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
site_sub <- comm_abund[facto$site_id=="105000010",]
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
site_sub <- comm_abund[facto$site_id=="105000011",]
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
site_sub <- comm_abund[facto$site_id=="105000012",]
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
site_sub <- comm_abund[facto$site_id=="105000013",]
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
site_sub <- comm_abund[facto$site_id=="105000014",]
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
site_sub <- comm_abund[facto$site_id=="105000015",]
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
site_sub <- comm_abund[facto$site_id=="105000016",]
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
site_sub <- comm_abund[facto$site_id=="105000017",]
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
site_sub <- comm_abund[facto$site_id=="105000018",]
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
site_sub <- comm_abund[facto$site_id=="105000019",]
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
site_sub <- comm_abund[facto$site_id=="105000020",]
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
site_sub <- comm_abund[facto$site_id=="106000001",]
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
site_sub <- comm_abund[facto$site_id=="106000002",]
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
site_sub <- comm_abund[facto$site_id=="106000003",]
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
site_sub <- comm_abund[facto$site_id=="106000004",]
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
site_sub <- comm_abund[facto$site_id=="106000005",]
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
site_sub <- comm_abund[facto$site_id=="106000006",]
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
site_sub <- comm_abund[facto$site_id=="106000007",]
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
site_sub <- comm_abund[facto$site_id=="106000008",]
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
site_sub <- comm_abund[facto$site_id=="106000009",]
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
site_sub <- comm_abund[facto$site_id=="110000001",]
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
site_sub <- comm_abund[facto$site_id=="110000002",]
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
site_sub <- comm_abund[facto$site_id=="110000003",]
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
site_sub <- comm_abund[facto$site_id=="110000004",]
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
site_sub <- comm_abund[facto$site_id=="110000005",]
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
site_sub <- comm_abund[facto$site_id=="110000006",]
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
site_sub <- comm_abund[facto$site_id=="110000007",]
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
site_sub <- comm_abund[facto$site_id=="110000008",]
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
site_sub <- comm_abund[facto$site_id=="110000009",]
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
site_sub <- comm_abund[facto$site_id=="110000010",]
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
site_sub <- comm_abund[facto$site_id=="110000011",]
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
site_sub <- comm_abund[facto$site_id=="110000012",]
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
site_sub <- comm_abund[facto$site_id=="110000013",]
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
site_sub <- comm_abund[facto$site_id=="110000014",]
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
site_sub <- comm_abund[facto$site_id=="110000015",]
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
site_sub <- comm_abund[facto$site_id=="110000016",]
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
site_sub <- comm_abund[facto$site_id=="112000001",]
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
site_sub <- comm_abund[facto$site_id=="113000001",]
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
site_sub <- comm_abund[facto$site_id=="113000002",]
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
site_sub <- comm_abund[facto$site_id=="113000003",]
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
site_sub <- comm_abund[facto$site_id=="113000004",]
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
site_sub <- comm_abund[facto$site_id=="113000005",]
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
site_sub <- comm_abund[facto$site_id=="115000001",]
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
site_sub <- comm_abund[facto$site_id=="115000002",]
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
site_sub <- comm_abund[facto$site_id=="116000004",]
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
site_sub <- comm_abund[facto$site_id=="116000007",]
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
site_sub <- comm_abund[facto$site_id=="118000004",]
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
site_sub <- comm_abund[facto$site_id=="119000001",]
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
site_sub <- comm_abund[facto$site_id=="119000002",]
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
site_sub <- comm_abund[facto$site_id=="119000003",]
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
site_sub <- comm_abund[facto$site_id=="119000004",]
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
site_sub <- comm_abund[facto$site_id=="119000005",]
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
site_sub <- comm_abund[facto$site_id=="119000006",]
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
site_sub <- comm_abund[facto$site_id=="119000007",]
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
site_sub <- comm_abund[facto$site_id=="119000008",]
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
site_sub <- comm_abund[facto$site_id=="119000009",]
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
site_sub <- comm_abund[facto$site_id=="119000010",]
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
site_sub <- comm_abund[facto$site_id=="120000001",]
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
site_sub <- comm_abund[facto$site_id=="120000002",]
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
site_sub <- comm_abund[facto$site_id=="120000003",]
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
site_sub <- comm_abund[facto$site_id=="101000001",]
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
site_sub <- comm_abund[facto$site_id=="101000002",]
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
site_sub <- comm_abund[facto$site_id=="101000003",]
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
site_sub <- comm_abund[facto$site_id=="101000004",]
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
site_sub <- comm_abund[facto$site_id=="101000005",]
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
site_sub <- comm_abund[facto$site_id=="101000006",]
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
site_sub <- comm_abund[facto$site_id=="101000007",]
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
site_sub <- comm_abund[facto$site_id=="101000008",]
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
site_sub <- comm_abund[facto$site_id=="101000009",]
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
site_sub <- comm_abund[facto$site_id=="101000010",]
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
write.csv(tur_o,"RareCountries_FD.csv")
##