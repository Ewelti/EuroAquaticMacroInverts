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
  }
  return(temp)
  
})

countryTrends <- do.call(rbind,countryTrends)
names(countryTrends)[which(names(countryTrends)=="siteID")] <- "site_id"
saveRDS(countryTrends,file="outputs/stanTrends_site_level_movingaverages.rds")

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

