#combine trends across all countries afterwards

trendsDir <- "C:/Users/db40fysa/Dropbox/Git/ellen_outputs"

trendsFiles <- list.files(trendsDir)[!grepl("txt",list.files(trendsDir))]

countryTrends <- lapply(trendsFiles,function(x){
   
  temp <- readRDS(paste(trendsDir,x,sep="/"))
  
  #add on response from file name
  temp$Response <- strsplit(x,"__")[[1]][2]

  return(temp)

})

countryTrends <- do.call(rbind,countryTrends)
saveRDS(countryTrends,file="outputs/stanTrends_site_level.rds")
