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

#all responses?
unique(countryTrends$Response)
#we are missing alien_Abund, alien_SppRich, E10, EPT_Abund, EPT_SppRich,
#insect_Abund, insect_SppRich

missing <- c("alien_Abund", "alien_SppRich", "E10", "EPT_Abund","EPT_SppRich",
             "insect_Abund", "insect_SppRich")

#see which ones have some NAs
temp <- subset(countryTrends, propNAs>0)
table(temp$Response)

missing2 <- names(table(temp$Response))

allmissing <- sort(c(missing,missing2))
