##Set working directory
#setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/")

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
siteData <- unique(d2[,c("site","Longitude_X","Latitude_Y")])
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
##color for positive values
pal_pos <- colorRampPalette(c('aliceblue','midnightblue')) #slategray1
##color for negative values
pal_neg <- colorRampPalette(c('red4','mistyrose'))

#spprich map
sr_sites <- resp[!is.na(resp$spp_richness),]
sr_sites <- sr_sites[order(sr_sites$spp_richness),]

#break positive and neg values
sr_neg <- sr_sites$spp_richness[(1:length(sr_sites$spp_richness))[sr_sites$spp_richness <= 0]]
sr_pos <- sr_sites$spp_richness[(1:length(sr_sites$spp_richness))[sr_sites$spp_richness > 0]]

#This adds a column of color values
sr_col_neg <- pal_neg(50)[as.numeric(cut(sr_neg,breaks = 50))]
unique(sr_col_neg)
sr_col_pos <- pal_pos(50)[as.numeric(cut(sr_pos,breaks = 50))]
unique(sr_col_pos)
sr_sites$sr_col <- c(sr_col_neg, sr_col_pos)

tiff(filename = "SppRich_map.tiff", width = 7, height = 6, units = 'in', res = 600, compression = 'lzw')

par(mar=c(0,0,0,0))

newmap <- getMap(resolution = "high")

col <- rep("grey70", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="aliceblue",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(34, 70),
  asp = 1,lwd=wv
)

sr_sites_asc <- sr_sites[order(-sr_sites$spp_richness),]
points(sr_sites$Longitude_X[sr_sites$spp_richness > 0],sr_sites$Latitude_Y[sr_sites$spp_richness > 0],pch = 20,col = alpha(sr_sites$sr_col[sr_sites$spp_richness > 0],0.6),cex=0.8)

points(sr_sites_asc$Longitude_X[sr_sites_asc$spp_richness <= 0],sr_sites_asc$Latitude_Y[sr_sites_asc$spp_richness <= 0],pch = 20,col = alpha(sr_sites_asc$sr_col[sr_sites_asc$spp_richness <= 0],0.6),cex=0.8)

lg <- round(seq(min(sr_sites$spp_richness), max(sr_sites$spp_richness), by=((max(sr_sites$spp_richness)-min(sr_sites$spp_richness))/7)),digits=1)

lg2 <- c(round(min(sr_sites$spp_richness),digits=1), "","",0, "","",round(max(sr_sites$spp_richness),digits=1))

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8.2,co_leng)
x2 <- rep(-6,co_leng)
#legend(-9, 72.3,title="",legend=lg,col =c(pal_neg(4),pal_pos(4)),lty=1,lwd=1,bty="n",cex=0.9)
legend(-9, 72.3,title="",legend=lg2,col =c('midnightblue','white','white','white','white','white','red4'),lty=1,lwd=1,bty="n",cex=1)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-10.2,72.7,legend=c("Species Richness"),bty='n')

dev.off()
#######################
##########################