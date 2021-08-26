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


### moving averages ####

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

#check we have all data
TaskID <- read.csv("outputs/MovingAverage_TaskIDs.csv",as.is=T)
TaskID$Index <- interaction(TaskID$country,TaskID$StartYear)
countryTrends$Index <- interaction(countryTrends$country,countryTrends$StartYear)

TaskID$Index[!TaskID$Index %in% countryTrends$Index]
# no :(

#create new TaskID with missing ones
TaskID <- subset(TaskID, !Index %in% countryTrends$Index)
TaskID$TaskID <- 1:nrow(TaskID)
write.table(TaskID,"outputs/MovingAverage_TaskIDs2.csv",sep=",",row.names=FALSE)

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
sens <- cbind(season, sens)

saveRDS(sens,file="sensitiv_seasonDiff.rds")

## seasonTrends
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
sens <- cbind(season, sens)

saveRDS(sens,file="sensitiv_seasonTrends.rds")

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
sens <- cbind(taxonres, sens)

saveRDS(sens,file="sensitiv_taxonresDiff.rds")

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
sens <- cbind(taxonres, sens)

saveRDS(sens,file="sensitiv_taxonresTrends.rds")