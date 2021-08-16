#combine trends across all countries for the same response afterwards

trendsDir <- "C:/Users/db40fysa/Dropbox/Git/ellen_outputs"
trendsFiles <- list.files(trendsDir)[!grepl("txt",list.files(trendsDir))]
countryTrends <- lapply(trendsFiles,function(x){
   readRDS(paste(trendsDir,x,sep="/"))
})
countryTrends <- do.call(rbind,countryTrends)
saveRDS(countryTrends,file="outputs/stanTrends_spp_richness_site-level.rds")
