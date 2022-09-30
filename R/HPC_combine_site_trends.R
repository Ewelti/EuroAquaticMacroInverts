### whole time series ####

#combine trends across all countries afterwards

trendsDir <- "C:/Users/db40fysa/Dropbox/Git/ellen_outputs"

#make sure it is just the rds trend files
trendsFiles <- list.files(trendsDir)[!grepl("txt",list.files(trendsDir))]
trendsFiles <- list.files(trendsDir)[grepl(".RDS",list.files(trendsDir))]

countryTrends <- lapply(trendsFiles,function(x){
   
  temp <- readRDS(paste(trendsDir,x,sep="/"))
  
  #add on response from file name
  temp$Response <- strsplit(as.character(x),"__")[[1]][2]
  return(temp)

})

countryTrends <- do.call(rbind,countryTrends)
names(countryTrends)[which(names(countryTrends)=="siteID")] <- "site_id"
saveRDS(countryTrends,file="outputs/stanTrends_site_level.rds")

#check we have all data

#get lists of tasks
TaskID <- read.csv("outputs/ResponseTrends_TaskIDs.csv",as.is=T)
sort(unique(TaskID$Response))
#yes!!

### updated models #####

#use beta for F_to
#use zero-inflated for AlienSppRich, EPTSppRich, EPTAbund - not use, since still doesnt work well when weighted
#the normal models for the rest

# standard model - for the rest
trendsDir <- "C:/Users/db40fysa/Dropbox/Git/ellen_outputs/siteLevel_updated"
trendsFiles <- allFiles <- list.files(trendsDir)[!grepl("txt",list.files(trendsDir))]

countryTrends <- lapply(trendsFiles,function(x){
  
  temp <- readRDS(paste(trendsDir,x,sep="/"))
  
  #add on response from file name
  temp$File <- x
  temp$Response <- strsplit(as.character(x),"__")[[1]][2]
  return(temp)
  
})
countryTrends <- do.call(rbind,countryTrends)

#beta models
trendsDir <- "C:/Users/db40fysa/Dropbox/Git/ellen_outputs/beta"
trendsFiles_Beta <- list.files(trendsDir)[!grepl("txt",list.files(trendsDir))]

countryTrends_Beta <- lapply(trendsFiles_Beta,function(x){
  
  temp <- readRDS(paste(trendsDir,x,sep="/"))
  
  #add on response from file name
  temp$File <- x
  temp$Response <- strsplit(as.character(x),"__")[[1]][2]
  return(temp)
  
})
countryTrends_Beta <- do.call(rbind,countryTrends_Beta)

#zero inflated models 
# trendsDir <- "C:/Users/db40fysa/Dropbox/Git/ellen_outputs/zero"
# trendsFiles_Zero <- list.files(trendsDir)[!grepl("txt",list.files(trendsDir))]
# 
# countryTrends_Zero <- lapply(trendsFiles_Zero,function(x){
# 
#   temp <- readRDS(paste(trendsDir,x,sep="/"))
# 
#   #add on response from file name
#   temp$File <- x
#   temp$Response <- strsplit(as.character(x),"__")[[1]][2]
#   return(temp)
# 
# })
# countryTrends_Zero <- do.call(rbind,countryTrends_Zero)


#combine all
#countryTrends <- subset(countryTrends, !File %in% countryTrends_Zero$File)
countryTrends <- subset(countryTrends, !File %in% countryTrends_Beta$File)
#countryTrends <- rbind(countryTrends,countryTrends_Beta,countryTrends_Zero)
countryTrends <- rbind(countryTrends,countryTrends_Beta)
names(countryTrends)[which(names(countryTrends)=="siteID")] <- "site_id"
saveRDS(countryTrends,file="outputs/stanTrends_site_level.rds")

##############################################
#### jackknife trends ####
trendsDir <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/meta_country_jacknife"
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/meta_country_jacknife")

ma <- readRDS("metaanalysis_country_jacknife_abund_nativeSpp.rds")
head(ma)

trendsFiles <- list.files(trendsDir)[grepl(".rds",list.files(trendsDir))]

countryTrends <- lapply(trendsFiles,function(x){
   
  temp <- readRDS(paste(trendsDir,x,sep="/"))
  
  #add on response from file name
  name <- gsub(".rds","", x)
  temp$Response <- gsub("metaanalysis_country_jacknife_","", name)
  return(temp)

})

countryTrends <- do.call(rbind,countryTrends)

# add country removed
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs")
siteData <- read.csv("All_siteLevel_and_glmOutput.csv", header=T) 
allCountries <- sort(unique(siteData$Country))
co <- rep(allCountries,20)
countryTrends$country <- co
head(countryTrends)
write.csv(countryTrends,file="JackknifeTrends.csv")

################################################

### site-level moving averages ####

trendsDir <- "C:/Users/db40fysa/Dropbox/Git/ellen_outputs/movingaverage"

#make sure it is just the rds trend files
trendsFiles <- list.files(trendsDir)[!grepl("txt",list.files(trendsDir))]
trendsFiles <- list.files(trendsDir)[grepl(".RDS",list.files(trendsDir))]

countryTrends <- lapply(trendsFiles,function(x){
  
  temp <- readRDS(paste(trendsDir,x,sep="/"))
  
  #add on response from file name
  if(nrow(temp)>0){
    temp$Response <- strsplit(as.character(x),"__")[[1]][2]
    temp$StartYear <- strsplit(as.character(x),"__")[[1]][4]
    temp$StartYear <- gsub(".RDS","",temp$StartYear)
    #temp$country <- strsplit(as.character(x),"__")[[1]][3]
  }
  return(temp)
  
})

countryTrends <- do.call(rbind,countryTrends)
names(countryTrends)[which(names(countryTrends)=="siteID")] <- "site_id"
saveRDS(countryTrends,file="outputs/stanTrends_site_level_movingaverages.rds")

### site-level MA higher threshold ####

TaskIDs <- read.csv("outputs/MovingAverageHigherThreshold_TaskIDs.csv")
originalTasks <- paste(TaskIDs$country, TaskIDs$StartYear, sep="__")
length(unique(originalTasks))#156 per response

modelFiles <- list.files("~/Dropbox/Collabs/Ellen/Ellen_HTMV")
modelTasks <- gsub("trendsHTMV__", "", modelFiles)
modelTasks <- gsub(".RDS", "", modelTasks)
modelResponse <- sapply(modelTasks, function(x) strsplit(x,"__")[[1]][1])
table(modelResponse)#all present!!

#combine all
countryTrends <- lapply(modelFiles,function(x){
  
  temp <- readRDS(paste("~/Dropbox/Collabs/Ellen/Ellen_HTMV",x,sep="/"))
  
  #add on response from file name
  if(nrow(temp)>0){
    temp$Response <- strsplit(as.character(x),"__")[[1]][2]
    temp$StartYear <- strsplit(as.character(x),"__")[[1]][4]
    temp$StartYear <- gsub(".RDS","",temp$StartYear)
    #temp$country <- strsplit(as.character(x),"__")[[1]][3]
  }
  return(temp)
  
})

countryTrends <- do.call(rbind,countryTrends)
names(countryTrends)[which(names(countryTrends)=="siteID")] <- "site_id"
saveRDS(countryTrends,file="outputs/stanTrends_site_level_movingaveragesHTMV.rds")

### moving average yr syntheses! #####
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/movingaverage_meta")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/movingaverage_meta"

require(data.table)
library(brms)

files = list.files(path = path, pattern = '\\.rds$')

dat_list = lapply(files, function(x){
  fit <- readRDS(x)
  fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
  fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
  fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
  fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
  fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
  nam <- gsub("metaanalysis_movingaverage_","", x)
  name <- gsub(".rds","", nam)
  name <- gsub("_1","__1", name)
  name <- gsub("_2","__2", name)
  response <- strsplit(as.character(name),"__")[[1]][1]
  year <- strsplit(as.character(name),"__")[[1]][2]
  fixed <- list(Response=response, StartYear=year, fixed_995[,1:4], fixed_975[,3:4],
                fixed_95[,3:4],fixed_90[,3:4])
  fixed <-data.frame(lapply(fixed, function(x) t(data.frame(x))))
  return(fixed)
})

MovAve <- do.call(rbind.data.frame, dat_list)
rownames(MovAve) <- NULL
head(MovAve)
##
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

ma <- readRDS("outputs/stanTrends_site_level_movingaverages.rds")
head(ma)
library(data.table)
unique(ma$Response)

#get site metadata
d1 <- read.csv("outputs/All_siteLevel_and_glmOutput.csv", header=T)
d1$site_id <- d1$site
head(d1)
siteData <- unique(d1[,c("site_id","study_id","Country","season","TaxonomicRes")])
ma <- merge(siteData,ma,by="site_id")

subs <- subset(ma, ma$Response == "abundance")
DT <- data.table(subs)
sitecount_ab <- DT[, .(site_num = length(unique(site_id))), by = StartYear]
sitecount_ab
COcount_ab <- DT[, .(country_num = length(unique(Country))), by = StartYear]
COcount_ab

subs_al <- subset(ma, ma$Response == "alien_SppRich")
DT_al <- data.table(subs_al)
sitecount_al <- DT_al[, .(site_num_Alien = length(unique(site_id))), by = StartYear]
sitecount_al

MoAv1 <- merge(MovAve,sitecount_ab,by="StartYear", all=T)
MoAv2 <- merge(MoAv1,sitecount_al,by="StartYear", all=T)
MoAv3 <- merge(MoAv2,COcount_ab,by="StartYear", all=T)
MoAv4 <- subset(MoAv3,StartYear >1989)
head(MoAv4)
##
write.csv(MoAv4, "outputs/movingAve_YrEsts.csv")

##############################################

### high threshold moving average yr syntheses! #####
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/HTMW1")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/HTMW1"

require(data.table)
library(brms)

files = list.files(path = path, pattern = '\\.rds$')

dat_list = lapply(files, function(x){
  fit <- readRDS(x)
  fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
  fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
  fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
  fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
  fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
  nam <- gsub("metaanalysis_movingaverage_higherthreshold","", x)
  name <- gsub(".rds","", nam)
  name <- gsub("_1","__1", name)
  name <- gsub("_2","__2", name)
  response <- strsplit(as.character(name),"__")[[1]][1]
  year <- strsplit(as.character(name),"__")[[1]][2]
  fixed <- list(Response=response, StartYear=year, fixed_995[,1:4], fixed_975[,3:4],
                fixed_95[,3:4],fixed_90[,3:4])
  fixed <-data.frame(lapply(fixed, function(x) t(data.frame(x))))
  return(fixed)
})
head(dat_list)

MovAve <- do.call(rbind.data.frame, dat_list)
head(MovAve)
unique(MovAve$Response)
rownames(MovAve) <- NULL
##
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

htma <- readRDS("outputs/stanTrends_site_level_movingaveragesHTMV.rds")
head(htma)
library(data.table)

#get site metadata
d1 <- read.csv("outputs/All_siteLevel_and_glmOutput.csv", header=T)
d1$site_id <- d1$site
head(d1)
siteData <- unique(d1[,c("site_id","study_id","Country","season","TaxonomicRes")])
htma <- merge(siteData,htma,by="site_id")

head(htma)

library(data.table)
subs <- subset(htma, htma$Response == "abundance")
DT <- data.table(subs)
sitecount_ab <- DT[, .(site_num = length(unique(site_id))), by = StartYear]
sitecount_ab
COcount_ab <- DT[, .(country_num = length(unique(Country))), by = StartYear]
COcount_ab

MoAv1 <- merge(MovAve,sitecount_ab,by="StartYear", all=T)
MoAv1 <- merge(MoAv1,COcount_ab,by="StartYear", all=T)
MoAv1 <- subset(MoAv1,StartYear >1999)
head(MoAv1)
##
write.csv(MoAv1, "outputs/HighThresholdMovingAve1_YrEsts.csv")
##

##############################################

### high threshold moving average yr syntheses2 #####
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/HTMW2")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/HTMW2"

require(data.table)
library(brms)

files = list.files(path = path, pattern = '\\.rds$')

dat_list = lapply(files, function(x){
  fit <- readRDS(x)
  fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
  fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
  fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
  fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
  fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
  nam <- gsub("metaanalysis_movingaverage_higherthreshold2_","", x)
  name <- gsub(".rds","", nam)
  name <- gsub("_1","__1", name)
  name <- gsub("_2","__2", name)
  response <- strsplit(as.character(name),"__")[[1]][1]
  year <- strsplit(as.character(name),"__")[[1]][2]
  fixed <- list(Response=response, StartYear=year, fixed_995[,1:4], fixed_975[,3:4],
                fixed_95[,3:4],fixed_90[,3:4])
  fixed <-data.frame(lapply(fixed, function(x) t(data.frame(x))))
  return(fixed)
})
head(dat_list)

MovAve <- do.call(rbind.data.frame, dat_list)
head(MovAve)
unique(MovAve$Response)
rownames(MovAve) <- NULL
##
#######################get number of sites per year
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")
library(lubridate)

##attach data
sites <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T) # change file name according to the time series to be analyzed
#delete missing data rows
DATA2 <- sites[!is.na(sites$site_id_wMissing),]

#count number of sites per year
sitecount <- aggregate(site_id ~ year, data = DATA2, FUN = length)

#subset by year
sites_later <- subset(DATA2,year > 1989)
head(sites_later)

#count number of sampling years per site
yearcount <- aggregate(year ~ site_id, data = sites_later, FUN = length)
head(yearcount)

#subset for sites with more years
siteslater_long <- yearcount[which(yearcount$year >= 20),]
head(siteslater_long)
nrow(siteslater_long)

###### get data to match with sites to count sites per year
response_stan <- readRDS("outputs/stanTrends_site_level_movingaverages.rds")

#subset estimates to sites selected for HTMW2
HTMW2sub <- subset(response_stan, site_id %in% siteslater_long$site_id)
head(HTMW2sub)
length(unique(HTMW2sub$site_id))

#get site metadata
d1 <- read.csv("outputs/All_siteLevel_and_glmOutput.csv", header=T)
d1$site_id <- d1$site
head(d1)
siteData <- unique(d1[,c("site_id","study_id","Country","season","TaxonomicRes")])
HTMW2sub <- merge(siteData,HTMW2sub,by="site_id")

head(HTMW2sub)
length(unique(HTMW2sub$site_id))
length(unique(HTMW2sub$Country))

library(data.table)
subs <- subset(HTMW2sub, HTMW2sub$Response == "abundance")
DT <- data.table(subs)
sitecount_ab <- DT[, .(site_num = length(unique(site_id))), by = StartYear]
sitecount_ab
COcount_ab <- DT[, .(country_num = length(unique(Country))), by = StartYear]
COcount_ab

MoAv1 <- merge(MovAve,sitecount_ab,by="StartYear", all=T)
MoAv1 <- merge(MoAv1,COcount_ab,by="StartYear", all=T)
MoAv1 <- subset(MoAv1,StartYear >1989)
head(MoAv1)

##
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")
write.csv(MoAv1, "outputs/HighThresholdMovingAve2_YrEsts.csv")
##
##############################################

#### Combine HT moving window trends
trendsDir <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/HTMW_siteLevel"

trendsFiles <- list.files(trendsDir)[grepl(".RDS",list.files(trendsDir))]

head(trendsFiles[1])

HTMWsiteTrends <- lapply(trendsFiles,function(x){
  
  temp <- data.frame(readRDS(paste(trendsDir,x,sep="/")))
  
  #add on response from file name
  temp$Response <- strsplit(as.character(x),"__")[[1]][2]
  temp$Country <- strsplit(as.character(x),"__")[[1]][3]
  temp$year <- strsplit(as.character(x),"__")[[1]][4]
  temp$year <- gsub(".RDS","",temp$year)
  return(temp)
  
})

HTMWsiteTrends <- do.call(rbind,HTMWsiteTrends)
head(HTMWsiteTrends)
saveRDS(HTMWsiteTrends,file="HTMWsiteTrends.rds")

#####################################################
##############################################

### species level sites only moving average #####
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/movingaverage_spplevelOnly")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/movingaverage_spplevelOnly"

require(data.table)
library(brms)

files = list.files(path = path, pattern = '\\.rds$')

dat_list = lapply(files, function(x){
  fit <- readRDS(x)
  fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
  fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
  fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
  fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
  fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
  nam <- gsub("metaanalysis_movingaverage_specieslevel_","", x)
  name <- gsub(".rds","", nam)
  name <- gsub("_1","__1", name)
  name <- gsub("_2","__2", name)
  response <- strsplit(as.character(name),"__")[[1]][1]
  year <- strsplit(as.character(name),"__")[[1]][2]
  fixed <- list(Response=response, StartYear=year, fixed_995[,1:4], fixed_975[,3:4],
                fixed_95[,3:4],fixed_90[,3:4])
  fixed <-data.frame(lapply(fixed, function(x) t(data.frame(x))))
  return(fixed)
})
head(dat_list)

MovAve <- do.call(rbind.data.frame, dat_list)
head(MovAve)
unique(MovAve$Response)
rownames(MovAve) <- NULL
##
#######################get number of sites per year
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")
library(lubridate)

##attach data
sites <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv", header=T) # change file name according to the time series to be analyzed
#delete missing data rows
DATA2 <- sites[!is.na(sites$site_id_wMissing),]
#subset for sites with spp level data
sppD <- DATA2[which(DATA2$TaxonomicRes == "species"),]
head(sppD)

#count number of sites per year
sitecount <- aggregate(site_id ~ year, data = sppD, FUN = length)
head(sitecount)

#count number of sampling years per site
yearcount <- aggregate(year ~ site_id, data = sppD, FUN = length)
head(yearcount)

###### get data to match with sites to count sites per year
response_stan <- readRDS("outputs/stanTrends_site_level_movingaverages.rds")

#subset estimates to sites selected for spp level only
SppSub <- subset(response_stan, site_id %in% sppD$site_id)

#get site metadata
d1 <- read.csv("outputs/All_siteLevel_and_glmOutput.csv", header=T)
d1$site_id <- d1$site
head(d1)
siteData <- unique(d1[,c("site_id","study_id","Country","season","TaxonomicRes")])
SppSub <- merge(siteData,SppSub,by="site_id")

head(SppSub)
length(unique(SppSub$site_id))
length(unique(SppSub$Country))

library(data.table)
subs <- subset(SppSub, SppSub$Response == "abundance")
DT <- data.table(subs)
sitecount_ab <- DT[, .(site_num = length(unique(site_id))), by = StartYear]
sitecount_ab
COcount_ab <- DT[, .(country_num = length(unique(Country))), by = StartYear]
COcount_ab

MoAv1 <- merge(MovAve,sitecount_ab,by="StartYear", all=T)
MoAv1 <- merge(MoAv1,COcount_ab,by="StartYear", all=T)
MoAv1 <- subset(MoAv1,StartYear >1989)
head(MoAv1)

##
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")
write.csv(MoAv1, "outputs/SppLevelMovAve_YrEsts.csv")
##
############################################################################################
##############################################

#### Combine trends for sensitivity analysis splitting the dataset by taxonomic resolution
sensDir <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/split_taxon_sensitivity"

sensFiles <- list.files(sensDir)[grepl(".rds",list.files(sensDir))]
head(sensFiles[1])

sensTrends <- lapply(sensFiles,function(x){
  
  temp <- data.frame(readRDS(paste(sensDir,x,sep="/")))
  
  #rename files to get rid of "_"
  x <- gsub("unweighted_","unweighted", x)
  x <- gsub("spp_richness","SppRich", x)
  x <- gsub("spp_rich_rare","SppRichRare", x)
  x <- gsub("SppRich_nativeSpp","SppRichNativeSpp", x)
  x <- gsub("abund_nativeSpp","abundNativeSpp", x)
  x <- gsub("alien_Abund","alienAbund", x)
  x <- gsub("alien_SppRich","alienSppRich", x)
  x <- gsub("EPT_Abund","EPTAbund", x)
  x <- gsub("EPT_SppRich","EPTSppRich", x)
  x <- gsub("F_to","Fto", x)
  x <- gsub("insect_Abund","insectAbund", x)
  x <- gsub("insect_SppRich","insectSppRich", x)

  #add on response from file name
  temp$resolution <- strsplit(as.character(x),"_")[[1]][2]
  temp$resolution <- gsub("taxonres","",temp$resolution)
  temp$response <- strsplit(as.character(x),"_")[[1]][3]
  temp$response <- gsub(".rds","",temp$response)
  return(temp)
  
})

sensTrends <- do.call(rbind,sensTrends)
head(sensTrends)
write.csv(sensTrends,file="SplitSensitivityTrends.csv")

##################################################################

### sensitivity analysis ####

library(raster)
require(data.table)

## seasonDiff in separate folder
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonDiff")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonDiff"

firstfile <- readRDS("fixef_seasonDiff_abundance.rds")
num_files <- length(list.files(path))
n.res <- 1:num_files
labs <- rownames(firstfile)
season <- rep(labs,length(n.res))

files = list.files(path = path, pattern = '\\.rds$')
dat_list = lapply(files, function (x) data.table(readRDS(x)))
trendsFiles <- gsub("fixef_seasonDiff_","", files)
trendsFiles <- gsub(".rds","", trendsFiles)
names(dat_list) <- trendsFiles
sens <- as.data.frame(dplyr::bind_rows(dat_list, .id = "Response"))
sens1 <- cbind(season, sens)

setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonDiff_nonNative")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonDiff_nonNative"

firstfile <- readRDS("fixef_seasonDiff_abund_nativeSpp.rds")
num_files <- length(list.files(path))
n.res <- 1:num_files
labs <- rownames(firstfile)
season <- rep(labs,length(n.res))

files = list.files(path = path, pattern = '\\.rds$')
dat_list = lapply(files, function (x) data.table(readRDS(x)))
trendsFiles <- gsub("fixef_seasonDiff_","", files)
trendsFiles <- gsub(".rds","", trendsFiles)
names(dat_list) <- trendsFiles
sens <- as.data.frame(dplyr::bind_rows(dat_list, .id = "Response"))
sens2 <- cbind(season, sens)

senss <- rbind(sens1, sens2)

saveRDS(senss,file="sensitiv_seasonDiff.rds")

### seasonTrends #####
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonTrends")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonTrends"

firstfile <- readRDS("fixef_seasonTrends_abundance.rds")
num_files <- length(list.files(path))
n.res <- 1:num_files
labs <- rownames(firstfile)
season <- rep(labs,length(n.res))

files = list.files(path = path, pattern = '\\.rds$')
dat_list = lapply(files, function (x) data.table(readRDS(x)))
trendsFiles <- gsub("fixef_seasonTrends_","", files)
trendsFiles <- gsub(".rds","", trendsFiles)
names(dat_list) <- trendsFiles
sens <- as.data.frame(dplyr::bind_rows(dat_list, .id = "Response"))
sens1 <- cbind(season, sens)

setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonTrends_nonNative")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonTrends_nonNative"

firstfile <- readRDS("fixef_seasonTrends_abund_nativeSpp.rds")
num_files <- length(list.files(path))
n.res <- 1:num_files
labs <- rownames(firstfile)
season <- rep(labs,length(n.res))

files = list.files(path = path, pattern = '\\.rds$')
dat_list = lapply(files, function (x) data.table(readRDS(x)))
trendsFiles <- gsub("fixef_seasonTrends_","", files)
trendsFiles <- gsub(".rds","", trendsFiles)
names(dat_list) <- trendsFiles
sens <- as.data.frame(dplyr::bind_rows(dat_list, .id = "Response"))
sens2 <- cbind(season, sens)

senss <- rbind(sens1, sens2)

saveRDS(senss,file="sensitiv_seasonTrends.rds")

## taxonresDiff
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresDiff")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresDiff"

firstfile <- readRDS("fixef_taxonresDiff_abundance.rds")
num_files <- length(list.files(path))
n.res <- 1:num_files
labs <- rownames(firstfile)
taxonres <- rep(labs,length(n.res))

files = list.files(path = path, pattern = '\\.rds$')
dat_list = lapply(files, function (x) data.table(readRDS(x)))
trendsFiles <- gsub("fixef_taxonresDiff_","", files)
trendsFiles <- gsub(".rds","", trendsFiles)
names(dat_list) <- trendsFiles
sens <- as.data.frame(dplyr::bind_rows(dat_list, .id = "Response"))
sens1 <- cbind(taxonres, sens)

setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresDiff_nonNative")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresDiff_nonNative"

firstfile <- readRDS("fixef_taxonresDiff_abund_nativeSpp.rds")
num_files <- length(list.files(path))
n.res <- 1:num_files
labs <- rownames(firstfile)
taxonres <- rep(labs,length(n.res))

files = list.files(path = path, pattern = '\\.rds$')
dat_list = lapply(files, function (x) data.table(readRDS(x)))
trendsFiles <- gsub("fixef_taxonresDiff_","", files)
trendsFiles <- gsub(".rds","", trendsFiles)
names(dat_list) <- trendsFiles
sens <- as.data.frame(dplyr::bind_rows(dat_list, .id = "Response"))
sens2 <- cbind(taxonres, sens)

senss <- rbind(sens1, sens2)

saveRDS(senss,file="sensitiv_taxonresDiff.rds")

## taxonresTrends
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresTrends")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresTrends"

firstfile <- readRDS("fixef_taxonresTrends_abundance.rds")
num_files <- length(list.files(path))
n.res <- 1:num_files
labs <- rownames(firstfile)
taxonres <- rep(labs,length(n.res))

files = list.files(path = path, pattern = '\\.rds$')
dat_list = lapply(files, function (x) data.table(readRDS(x)))
trendsFiles <- gsub("fixef_taxonresTrends_","", files)
trendsFiles <- gsub(".rds","", trendsFiles)
names(dat_list) <- trendsFiles
sens <- as.data.frame(dplyr::bind_rows(dat_list, .id = "Response"))
sens1 <- cbind(taxonres, sens)

setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresTrends_nonNative")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresTrends_nonNative"

firstfile <- readRDS("fixef_taxonresTrends_abund_nativeSpp.rds")
num_files <- length(list.files(path))
n.res <- 1:num_files
labs <- rownames(firstfile)
taxonres <- rep(labs,length(n.res))

files = list.files(path = path, pattern = '\\.rds$')
dat_list = lapply(files, function (x) data.table(readRDS(x)))
trendsFiles <- gsub("fixef_taxonresDiff_","", files)
trendsFiles <- gsub(".rds","", trendsFiles)
names(dat_list) <- trendsFiles
sens <- as.data.frame(dplyr::bind_rows(dat_list, .id = "Response"))
sens2 <- cbind(taxonres, sens)

senss <- rbind(sens1, sens2)

saveRDS(senss,file="sensitiv_taxonresTrends.rds")

##########################################
#### one stage decadal trends #####
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/onestage_decadeEffect")
path <- "C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/onestage_decadeEffect"

htma <- readRDS("onestage_decadeEffect_abund_nativeSpp.rds")
head(htma)

require(data.table)
library(brms)

files = list.files(path = path, pattern = '\\.rds$')

DecadeTrends <- lapply(files,function(x){
  
  temp <- data.frame(readRDS(paste(path,x,sep="/")))
  
  #add on response from file name
  name <- gsub(".rds","", x)
  temp$Response <- strsplit(as.character(name),"onestage_decadeEffect_")[[1]][2]
  temp$var <- rownames(temp) 
  return(temp)
  
})

DecadeTrends <- do.call(rbind,DecadeTrends)
rownames(DecadeTrends) <- NULL
DecadeTrends <- DecadeTrends[, c(6, 5, 1, 2, 3, 4)]
head(DecadeTrends)
write.csv(DecadeTrends,file="DecadeTrends.csv")
###########################################################

### country level trends ####

trendsDir <- "C:/Users/db40fysa/Dropbox/Git/ellen_outputs/metacountries"

#make sure it is just the rds trend files
trendsFiles <- list.files(trendsDir)[!grepl("txt",list.files(trendsDir))]
trendsFiles <- list.files(trendsDir)[grepl(".rds",list.files(trendsDir))]

countryTrends <- lapply(trendsFiles,function(x){
  
  temp <- data.frame(readRDS(paste(trendsDir,x,sep="/")))
  
  #add on response from file name
  temp$Response <- strsplit(as.character(x),"__")[[1]][2]
  temp$Country <- strsplit(as.character(x),"__")[[1]][3]
  temp$Country <- gsub(".rds","",temp$Country)
  return(temp)
  
})

countryTrends <- do.call(rbind,countryTrends)
saveRDS(countryTrends,file="outputs/stanTrends_country_level.rds")

library(ggplot2)
library(ggthemes)

ggplot(countryTrends)+
  geom_pointrange(aes(x=Country, y=Estimate, ymin=Q2.5, ymax=Q97.5))+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  facet_wrap(~Response)+
  theme_few()

#####################################################

### moving average yr syntheses from splits by latitude! #####
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/movingaverage_meta_split")
path <- "C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/movingaverage_meta_split"

require(data.table)
library(brms)

files = list.files(path = path, pattern = '\\.rds$')

dat_list = lapply(files, function(x){
  fit <- readRDS(x)
  fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
  fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
  fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
  fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
  fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
  nam <- gsub("metaanalysis_movingaverage_","", x)
  name <- gsub(".rds","", nam)
  name <- gsub("_1","__1", name)
  name <- gsub("_2","__2", name)
  response <- strsplit(as.character(name),"__")[[1]][1]
  year <- strsplit(as.character(name),"__")[[1]][2]
  fixed <- list(Response=response, StartYear=year, fixed_995[,1:4], fixed_975[,3:4],
                fixed_95[,3:4],fixed_90[,3:4])
  fixed <-data.frame(lapply(fixed, function(x) t(data.frame(x))))
  return(fixed)
})

MovAve <- do.call(rbind.data.frame, dat_list)
head(MovAve)
##
write.csv(MovAve, "outputs/movingAve_split_YrEsts.csv")
###############################################################
