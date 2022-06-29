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

#some missing - check

### moving average yr syntheses! #####
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/movingaverage_meta")
path <- "C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/movingaverage_meta"

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
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/")

ma <- readRDS("outputs/stanTrends_site_level_movingaverages.rds")
head(ma)
library(data.table)

subs <- subset(ma, ma$Response == "abundance")
DT <- data.table(subs)
sitecount_ab <- DT[, .(site_num_Abund = length(unique(site_id))), by = StartYear]
sitecount_ab

subs_sr <- subset(ma, ma$Response == "spp_richness")
DT_sr <- data.table(subs_sr)
sitecount_sr <- DT_sr[, .(site_num_SRich = length(unique(site_id))), by = StartYear]
sitecount_sr

subs_nat <- subset(ma, ma$Response == "SppRich_nativeSpp")
DT_nat <- data.table(subs_nat)
sitecount_nat <- DT_nat[, .(site_num_Native = length(unique(site_id))), by = StartYear]
sitecount_nat

subs_al <- subset(ma, ma$Response == "alien_SppRich")
DT_al <- data.table(subs_al)
sitecount_al <- DT_al[, .(site_num_Alien = length(unique(site_id))), by = StartYear]
sitecount_al

subs_FRic <- subset(ma, ma$Response == "FRic")
DT_FRic <- data.table(subs_FRic)
sitecount_FRic <- DT_FRic[, .(site_num_FRic = length(unique(site_id))), by = StartYear]
sitecount_FRic

subs_FRed <- subset(ma, ma$Response == "FRed")
DT_FRed <- data.table(subs_FRed)
sitecount_FRed <- DT_FRed[, .(site_num_FRed = length(unique(site_id))), by = StartYear]
sitecount_FRed

MoAv1 <- merge(MovAve,sitecount_ab,by="StartYear", all=T)
MoAv2 <- merge(MoAv1,sitecount_sr,by="StartYear", all=T)
MoAv3 <- merge(MoAv2,sitecount_nat,by="StartYear", all=T)
MoAv4 <- merge(MoAv3,sitecount_al,by="StartYear", all=T)
MoAv5 <- merge(MoAv4,sitecount_FRic,by="StartYear", all=T)
MoAv6 <- merge(MoAv5,sitecount_FRed,by="StartYear", all=T)
head(MoAv6)
##
write.csv(MoAv6, "outputs/movingAve_YrEsts.csv")

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

### sensitivity analysis ####

library(raster)
require(data.table)

## seasonDiff in separate folder
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonDiff")
path <- "C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonDiff"

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

setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonDiff_nonNative")
path <- "C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonDiff_nonNative"

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
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonTrends")
path <- "C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonTrends"

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

setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonTrends_nonNative")
path <- "C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/seasonTrends_nonNative"

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
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresDiff")
path <- "C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresDiff"

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

setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresDiff_nonNative")
path <- "C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresDiff_nonNative"

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
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresTrends")
path <- "C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresTrends"

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

setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresTrends_nonNative")
path <- "C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/outputs/Sensitivity/taxonresTrends_nonNative"

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
