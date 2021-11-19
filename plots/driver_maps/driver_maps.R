##Set working directory
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/")

# attach data
response_stan <- readRDS("outputs/stanTrends_site_level.rds")
head(response_stan)
unique(response_stan$Response)

library(tidyverse)
#pivot responses
response_stan_pivot <- response_stan %>%
                    select(c(Response,estimate,site_id)) %>%
                    pivot_wider(names_from = "Response",
                                values_from = "estimate")

#get site metadata
d2 <- read.csv("outputs/All_siteLevel_and_glmOutput.csv", header=T)
siteData <- unique(d2[,c("site","Longitude_X","Latitude_Y","ppt_Est", "tmax_Est", 
                         "ppt_mm_12moPrior", "tmax_C_12moPrior", 
                         "strahler_streamOrder","accumulation_atPoint", 
                         "elevation_atPoint","slope_mean", 
                         "N_Est", "N_mean", 
                         "urban_meanPerc_upstr", "crop_meanPerc_upstr",
                         "crop_Est", "urban_Est", 
                         "dam_impact_score_lessthan100km")])
colnames(siteData)[1] <- "site_id"
resp <- merge(siteData,response_stan_pivot,by="site_id")
head(resp)
#################################################################################

library(rworldxtra)
library(RColorBrewer)
library(maps)
library(rworldmap)
library(scales)

#####################
###################
getMap()$NAME

ddf = read.table(text="
country
'France' 
'Estonia'
'Cyprus'
'N. Cyprus'
'Spain'
'Sweden'
'Luxembourg'
'Bulgaria'
'Denmark'
'Germany'
'United Kingdom'
'Ireland'
'Norway'
'Czech Rep.'
'Italy'
'Hungary'
'Portugal'
'Austria'
'Netherlands'
'Switzerland'
'Finland'
'Latvia'
'Belgium'", header=TRUE)

######################
####################
#Create a function to generate a continuous color palette
pal_ <- colorRampPalette(c("darkred",'aliceblue','midnightblue')) #slategray1

pal_ <- colorRampPalette(c('darkred','gold', 'seagreen2','dodgerblue4','darkviolet'))
################### DAM map #############################
tiff(filename = "plots/driver_maps/Dam_map.tiff", width = 6.8, height = 5.8, units = 'in', res = 600, compression = 'lzw')

pal_ <- colorRampPalette(c('gold', 'seagreen2','dodgerblue4','darkviolet','darkred'))
resp$response <- resp$dam_impact_score_lessthan100km
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#This adds a column of color values
sr_sites$sr_col <- pal_(50)[as.numeric(cut(sr_sites$response,breaks = 50))]
unique(sr_sites$sr_col)

par(mar=c(0,0,0,0))

newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))#"white"
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70", #bg="white",border=col2,
  xlim = c(-10, 34),
  ylim = c(34, 70),
  asp = 1,lwd=wv
)

sr_sites_asc <- sr_sites[order(sr_sites$response),]
points(sr_sites_asc$Longitude_X,sr_sites_asc$Latitude_Y,pch = 20,col = alpha(sr_sites_asc$sr_col,0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=1)

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9, 72.3,title="",legend=lg,col =pal_(8),lty=1,lwd=1,bty="n",cex=0.85)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3.5)
legend(-10.2,72.7,legend=c("Dam Impact"),bty='n')

dev.off()
#######################
##########################


########################  CLIMATE  #####################
tiff(filename = "plots/driver_maps/Climate_maps.tiff", width = 11, height = 11, units = 'in', res = 600, compression = 'lzw')

par(mar=c(0,0,0,0),mfrow=c(2,2))
###################ppt_mm_12moPrior map#############################
pal_ <- colorRampPalette(c('darkred','gold', 'seagreen2','dodgerblue4','darkviolet'))
resp$response <- resp$ppt_mm_12moPrior
resp$logresponse <- log10(resp$response)
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#This adds a column of color values
sr_sites$sr_col <- pal_(50)[as.numeric(cut(sr_sites$logresponse,breaks = 50))]
sr_sites$log_sr_col <- pal_(50)[as.numeric(cut(sr_sites$logresponse,breaks = 50))]

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(30, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(sr_sites$response),]
points(sr_sites_asc$Longitude_X,sr_sites_asc$Latitude_Y,pch = 20,col = alpha(sr_sites_asc$log_sr_col,0.6),cex=1.5)

lg_log <- round(seq(min(sr_sites$logresponse), max(sr_sites$logresponse), by=((max(sr_sites$logresponse)-min(sr_sites$logresponse))/7)),digits=6)

##inverse log!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
lg <- round(10^lg_log,digits=0)

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9, 72.3,title="",legend=lg,col =pal_(8),lty=1,lwd=1,bty="n",cex=0.85)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-10.2,72.7,legend=c("Mean Precip (mm/yr)"),bty='n')

###################ppt_Est map#############################
##color for positive values
pal_pos <- colorRampPalette(c('aliceblue','midnightblue')) #slategray1
##color for negative values
pal_neg <- colorRampPalette(c('red4','mistyrose'))

resp$response <- resp$ppt_Est
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#break positive and neg values
sr_neg <- sr_sites$response[(1:length(sr_sites$response))[sr_sites$response <= 0]]
sr_pos <- sr_sites$response[(1:length(sr_sites$response))[sr_sites$response > 0]]

#This adds a column of color values
sr_col_neg <- pal_neg(50)[as.numeric(cut(sr_neg,breaks = 50))]
unique(sr_col_neg)
sr_col_pos <- pal_pos(50)[as.numeric(cut(sr_pos,breaks = 50))]
unique(sr_col_pos)
sr_sites$sr_col <- c(sr_col_neg, sr_col_pos)

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(30, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(-sr_sites$response),]
points(sr_sites$Longitude_X[sr_sites$response > 0],sr_sites$Latitude_Y[sr_sites$response > 0],pch = 20,col = alpha(sr_sites$sr_col[sr_sites$response > 0],0.6),cex=1.5)

points(sr_sites_asc$Longitude_X[sr_sites_asc$response <= 0],sr_sites_asc$Latitude_Y[sr_sites_asc$response <= 0],pch = 20,col = alpha(sr_sites_asc$sr_col[sr_sites_asc$response <= 0],0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=1)

lg2 <- c(round(min(sr_sites$response),digits=1), "","",0, "","",round(max(sr_sites$response),digits=1))

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8.2,co_leng)
x2 <- rep(-6,co_leng)
#legend(-9, 72.3,title="",legend=lg,col =c(pal_neg(4),pal_pos(4)),lty=1,lwd=1,bty="n",cex=0.94)
legend(-9, 72.3,title="",legend=lg2,col =c('midnightblue','white','white','white','white','white','red4'),lty=1,lwd=1,bty="n",cex=0.94)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-10.2,72.7,legend=c("Slope of Precip (mm/yr)"),bty='n')

###################tmax_C_12moPrior map#############################
pal_ <- colorRampPalette(c('dodgerblue4','seagreen2','gold','darkred','darkviolet'))
resp$response <- resp$tmax_C_12moPrior
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#This adds a column of color values
sr_sites$sr_col <- pal_(50)[as.numeric(cut(sr_sites$response,breaks = 50))]

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(30, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(sr_sites$response),]
points(sr_sites_asc$Longitude_X,sr_sites_asc$Latitude_Y,pch = 20,col = alpha(sr_sites_asc$sr_col,0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=1)

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9, 72.3,title="",legend=lg,col =pal_(8),lty=1,lwd=1,bty="n",cex=0.85)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-10.2,72.7,legend=c("Mean Tmax (°C/yr)"),bty='n')

###################tmax_Est map#############################
##color for positive values
pal_pos <- colorRampPalette(c('mistyrose','red4'))
##color for negative values
pal_neg <- colorRampPalette(c('midnightblue','aliceblue'))

resp$response <- resp$tmax_Est
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#break positive and neg values
sr_neg <- sr_sites$response[(1:length(sr_sites$response))[sr_sites$response <= 0]]
sr_pos <- sr_sites$response[(1:length(sr_sites$response))[sr_sites$response > 0]]

#This adds a column of color values
sr_col_neg <- pal_neg(50)[as.numeric(cut(sr_neg,breaks = 50))]
unique(sr_col_neg)
sr_col_pos <- pal_pos(50)[as.numeric(cut(sr_pos,breaks = 50))]
unique(sr_col_pos)
sr_sites$sr_col <- c(sr_col_neg, sr_col_pos)

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(30, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(-sr_sites$response),]
points(sr_sites$Longitude_X[sr_sites$response > 0],sr_sites$Latitude_Y[sr_sites$response > 0],pch = 20,col = alpha(sr_sites$sr_col[sr_sites$response > 0],0.6),cex=1.5)

points(sr_sites_asc$Longitude_X[sr_sites_asc$response <= 0],sr_sites_asc$Latitude_Y[sr_sites_asc$response <= 0],pch = 20,col = alpha(sr_sites_asc$sr_col[sr_sites_asc$response <= 0],0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=1)

lg2 <- c(round(min(sr_sites$response),digits=1), "",0,"","","",round(max(sr_sites$response),digits=1))

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8.2,co_leng)
x2 <- rep(-6,co_leng)
#legend(-9, 72.3,title="",legend=lg,col =c(pal_neg(4),pal_pos(4)),lty=1,lwd=1,bty="n",cex=0.9)
legend(-9, 72.3,title="",legend=lg2,col =c('midnightblue','white','white','white','white','white','red4'),lty=1,lwd=1,bty="n",cex=0.94)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-10.2,72.7,legend=c("Slope of Tmax (°C/yr)"),bty='n')

dev.off()
#######################
##########################


########################  LAND COVER  #####################
tiff(filename = "plots/driver_maps/Landcover_maps.tiff", width = 11, height = 6, units = 'in', res = 600, compression = 'lzw')

par(mar=c(0,0,0,0),mfrow=c(1,2))
pal_ <- colorRampPalette(c('seagreen2','gold','darkred'))#,'darkviolet'))
###################crop_meanPerc_upstr map#############################
resp$response <- resp$crop_meanPerc_upstr
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#This adds a column of color values
sr_sites$sr_col <- pal_(50)[as.numeric(cut(sr_sites$response,breaks = 50))]

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(30, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(sr_sites$response),]
points(sr_sites_asc$Longitude_X,sr_sites_asc$Latitude_Y,pch = 20,col = alpha(sr_sites_asc$sr_col,0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=1)
lg <- lg*100

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9, 72.3,title="",legend=lg,col =pal_(8),lty=1,lwd=1,bty="n",cex=0.7)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-10.2,72.7,legend=c("% Crop"),bty='n')
###############################################################
###################urban_meanPerc_upstr map#############################

resp$response <- resp$urban_meanPerc_upstr
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#This adds a column of color values
sr_sites$sr_col <- pal_(50)[as.numeric(cut(sr_sites$response,breaks = 50))]

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(30, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(sr_sites$response),]
points(sr_sites_asc$Longitude_X,sr_sites_asc$Latitude_Y,pch = 20,col = alpha(sr_sites_asc$sr_col,0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=1)
lg <- lg*100

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9, 72.3,title="",legend=lg,col =pal_(8),lty=1,lwd=1,bty="n",cex=0.7)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-10.2,72.7,legend=c("% Urban"),bty='n')

dev.off()
#######################


#######################
##########################


########################  STREAM CHARACTERISTICS  #####################
tiff(filename = "plots/driver_maps/StreamCharaterictics_maps.tiff", width = 11, height = 11, units = 'in', res = 600, compression = 'lzw')

par(mar=c(0,0,0,0),mfrow=c(2,2))
pal_ <- colorRampPalette(c('gold', 'seagreen2','dodgerblue4','darkviolet','darkred'))
###################strahler_streamOrder map#############################
resp$response <- resp$strahler_streamOrder
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#This adds a column of color values
sr_sites$sr_col <- pal_(50)[as.numeric(cut(sr_sites$response,breaks = 50))]

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(30, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(sr_sites$response),]
points(sr_sites_asc$Longitude_X,sr_sites_asc$Latitude_Y,pch = 20,col = alpha(sr_sites_asc$sr_col,0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=0)

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9, 72.3,title="",legend=lg,col =pal_(8),lty=1,lwd=1,bty="n",cex=0.85)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=11)
legend(-10.2,72.7,legend=c("Stream Order"),bty='n')

###################accumulation_atPoint map#############################
resp$response <- resp$accumulation_atPoint
resp$logresponse <- log10(resp$response)
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#This adds a column of color values
sr_sites$sr_col <- pal_(50)[as.numeric(cut(sr_sites$response,breaks = 50))]
sr_sites$log_sr_col <- pal_(50)[as.numeric(cut(sr_sites$logresponse,breaks = 50))]

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(30, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(sr_sites$response),]
points(sr_sites_asc$Longitude_X,sr_sites_asc$Latitude_Y,pch = 20,col = alpha(sr_sites_asc$log_sr_col,0.6),cex=1.5)

lg_log <- round(seq(min(sr_sites$logresponse), max(sr_sites$logresponse), by=((max(sr_sites$logresponse)-min(sr_sites$logresponse))/7)),digits=9)

##inverse log!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
lg <- round(10^lg_log,digits=2)

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9, 72.3,title="",legend=lg,col =pal_(8),lty=1,lwd=1,bty="n",cex=0.85)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=5)
legend(-10.2,72.7,legend=expression(paste("Flow Accumulation (",km^2,")", sep = "")),bty='n')

###################elevation_atPoint map#############################
resp$response <- resp$elevation_atPoint
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#This adds a column of color values
sr_sites$sr_col <- pal_(50)[as.numeric(cut(sr_sites$response,breaks = 50))]

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(30, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(sr_sites$response),]
points(sr_sites_asc$Longitude_X,sr_sites_asc$Latitude_Y,pch = 20,col = alpha(sr_sites_asc$sr_col,0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=0)

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9, 72.3,title="",legend=lg,col =pal_(8),lty=1,lwd=1,bty="n",cex=0.85)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-10.2,72.7,legend=c("Elevation (m)"),bty='n')

###################slope_mean map#############################
resp$response <- resp$slope_mean
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#This adds a column of color values
sr_sites$sr_col <- pal_(50)[as.numeric(cut(sr_sites$response,breaks = 50))]

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(30, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(sr_sites$response),]
points(sr_sites_asc$Longitude_X,sr_sites_asc$Latitude_Y,pch = 20,col = alpha(sr_sites_asc$sr_col,0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=1)
lg <- lg*100

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9, 72.3,title="",legend=lg,col =pal_(8),lty=1,lwd=1,bty="n",cex=0.85)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=4)
legend(-10.2,72.7,legend=c("Slope (%)"),bty='n')

dev.off()
#######################

#######################
##########################


########################  INORGANIC N  #####################
tiff(filename = "InorganicN_maps.tiff", width = 11, height = 6, units = 'in', res = 600, compression = 'lzw')

par(mar=c(0,0,0,0),mfrow=c(1,2))
###################N_mean map#############################
resp$response <- resp$N_mean
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#This adds a column of color values
sr_sites$sr_col <- pal_(50)[as.numeric(cut(sr_sites$response,breaks = 50))]

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(34, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(sr_sites$response),]
points(sr_sites_asc$Longitude_X,sr_sites_asc$Latitude_Y,pch = 20,col = alpha(sr_sites_asc$sr_col,0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=0)

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9, 72.3,title="",legend=lg,col =pal_(8),lty=1,lwd=1,bty="n",cex=0.7)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-10.2,72.7,legend=c("Mean Inorg. N (kg/ha/Spr. month)"),bty='n')

###################N_Est map#############################
##color for positive values
pal_pos <- colorRampPalette(c('mistyrose','red4'))
##color for negative values
pal_neg <- colorRampPalette(c('midnightblue','aliceblue'))

resp$response <- resp$N_Est
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#break positive and neg values
sr_neg <- sr_sites$response[(1:length(sr_sites$response))[sr_sites$response <= 0]]
sr_pos <- sr_sites$response[(1:length(sr_sites$response))[sr_sites$response > 0]]

#This adds a column of color values
sr_col_neg <- pal_neg(50)[as.numeric(cut(sr_neg,breaks = 50))]
unique(sr_col_neg)
sr_col_pos <- pal_pos(50)[as.numeric(cut(sr_pos,breaks = 50))]
unique(sr_col_pos)
sr_sites$sr_col <- c(sr_col_neg, sr_col_pos)

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(34, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(-sr_sites$response),]
points(sr_sites_asc$Longitude_X[sr_sites_asc$response <= 0],sr_sites_asc$Latitude_Y[sr_sites_asc$response <= 0],pch = 20,col = alpha(sr_sites_asc$sr_col[sr_sites_asc$response <= 0],0.6),cex=1.5)

points(sr_sites$Longitude_X[sr_sites$response > 0],sr_sites$Latitude_Y[sr_sites$response > 0],pch = 20,col = alpha(sr_sites$sr_col[sr_sites$response > 0],0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=1)

lg2 <- c(round(min(sr_sites$response),digits=1), "","","","","","0",round(max(sr_sites$response),digits=1))

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8.2,co_leng)
x2 <- rep(-6,co_leng)
#legend(-9, 72.3,title="",legend=lg,col =c(pal_neg(4),pal_pos(4)),lty=1,lwd=1,bty="n",cex=0.94)
legend(-9, 72.3,title="",legend=lg2,col =c('midnightblue','white','white','white','white','white','white','red4'),lty=1,lwd=1,bty="n",cex=0.7)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=4.5)
legend(-10.2,72.7,legend=c("Slope of Inorg. N (kg/ha/Spr. month)"),bty='n')


dev.off()
##
##ää





###################urban_Est map#############################
##color for positive values
pal_pos <- colorRampPalette(c('mistyrose','red4'))
##color for negative values
pal_neg <- colorRampPalette(c('midnightblue','aliceblue'))

resp$response <- resp$urban_Est
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#break positive and neg values
sr_neg <- sr_sites$response[(1:length(sr_sites$response))[sr_sites$response <= 0]]
sr_pos <- sr_sites$response[(1:length(sr_sites$response))[sr_sites$response > 0]]

#This adds a column of color values
sr_col_neg <- pal_neg(50)[as.numeric(cut(sr_neg,breaks = 50))]
unique(sr_col_neg)
sr_col_pos <- pal_pos(50)[as.numeric(cut(sr_pos,breaks = 50))]
unique(sr_col_pos)
sr_sites$sr_col <- c(sr_col_neg, sr_col_pos)

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(34, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(-sr_sites$response),]
points(sr_sites_asc$Longitude_X[sr_sites_asc$response <= 0],sr_sites_asc$Latitude_Y[sr_sites_asc$response <= 0],pch = 20,col = alpha(sr_sites_asc$sr_col[sr_sites_asc$response <= 0],0.6),cex=1.5)

points(sr_sites$Longitude_X[sr_sites$response > 0],sr_sites$Latitude_Y[sr_sites$response > 0],pch = 20,col = alpha(sr_sites$sr_col[sr_sites$response > 0],0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=1)

lg2 <- c(round(min(sr_sites$response),digits=1), "","",0, "","",round(max(sr_sites$response),digits=1))

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8.2,co_leng)
x2 <- rep(-6,co_leng)
#legend(-9, 72.3,title="",legend=lg,col =c(pal_neg(4),pal_pos(4)),lty=1,lwd=1,bty="n",cex=0.94)
legend(-9, 72.3,title="",legend=lg2,col =c('midnightblue','white','white','white','white','white','red4'),lty=1,lwd=1,bty="n",cex=0.94)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=4.5)
legend(-10.2,72.7,legend=c("Slope of % Urban"),bty='n')

##########################################################
###################################
###################crop_Est map#############################
##color for positive values
pal_pos <- colorRampPalette(c('mistyrose','red4'))
##color for negative values
pal_neg <- colorRampPalette(c('midnightblue','aliceblue'))

resp$response <- resp$tmax_Est
sr_sites <- resp[!is.na(resp$response),]
sr_sites <- sr_sites[order(sr_sites$response),]

#break positive and neg values
sr_neg <- sr_sites$response[(1:length(sr_sites$response))[sr_sites$response <= 0]]
sr_pos <- sr_sites$response[(1:length(sr_sites$response))[sr_sites$response > 0]]

#This adds a column of color values
sr_col_neg <- pal_neg(50)[as.numeric(cut(sr_neg,breaks = 50))]
unique(sr_col_neg)
sr_col_pos <- pal_pos(50)[as.numeric(cut(sr_pos,breaks = 50))]
unique(sr_col_pos)
sr_sites$sr_col <- c(sr_col_neg, sr_col_pos)

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

col <- rep("white", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(34, 70),
  asp = 1,lwd=wv)

sr_sites_asc <- sr_sites[order(-sr_sites$response),]
points(sr_sites$Longitude_X[sr_sites$response > 0],sr_sites$Latitude_Y[sr_sites$response > 0],pch = 20,col = alpha(sr_sites$sr_col[sr_sites$response > 0],0.6),cex=1.5)

points(sr_sites_asc$Longitude_X[sr_sites_asc$response <= 0],sr_sites_asc$Latitude_Y[sr_sites_asc$response <= 0],pch = 20,col = alpha(sr_sites_asc$sr_col[sr_sites_asc$response <= 0],0.6),cex=1.5)

lg <- round(seq(min(sr_sites$response), max(sr_sites$response), by=((max(sr_sites$response)-min(sr_sites$response))/7)),digits=1)

lg2 <- c(round(min(sr_sites$response),digits=1), "","",0, "","",round(max(sr_sites$response),digits=1))

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8.2,co_leng)
x2 <- rep(-6,co_leng)
#legend(-9, 72.3,title="",legend=lg,col =c(pal_neg(4),pal_pos(4)),lty=1,lwd=1,bty="n",cex=0.9)
legend(-9, 72.3,title="",legend=lg2,col =c('midnightblue','white','white','white','white','white','red4'),lty=1,lwd=1,bty="n",cex=0.94)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-10.2,72.7,legend=c("Slope of % Crop"),bty='n')


