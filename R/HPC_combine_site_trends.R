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
