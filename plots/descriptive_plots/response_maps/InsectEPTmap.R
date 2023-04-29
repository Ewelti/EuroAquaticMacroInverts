##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")

# attach data
response_stan <- readRDS("outputs/outputs_metaAnalysis/stanTrends_site_level.rds")
head(response_stan)
unique(response_stan$Response)

library(tidyverse)
#pivot responses
response_stan_pivot <- response_stan %>%
                    select(c(Response,estimate,site_id)) %>%
                    pivot_wider(names_from = "Response",
                                values_from = "estimate")

#get site metadata
d2 <- read.csv("outputs/All_siteLevel.csv", header=T)
siteData <- unique(d2[,c("site","Longitude_X","Latitude_Y")])
colnames(siteData)[1] <- "site_id"
resp <- merge(siteData,response_stan_pivot,by="site_id")
head(resp)

metricdata <- read.csv("outputs/All_indices_benthicMacroInverts_AllYears_alienzeros.csv")
head(metricdata)
#################################################################################
##install.packages("rworldmap")
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

tiff(filename = "plots/descriptive_plots/response_maps/InsectEPT_map.tiff", width = 11, height = 11, units = 'in', res = 600, compression = 'lzw')

par(mar=c(0,0,0,0),mfrow=c(2,2))
############################EPT_SppRich#######################################
#Create a function to generate a continuous color palette
##color for positive values
pal_pos <- colorRampPalette(c('aliceblue','midnightblue')) #slategray1
##color for negative values
pal_neg <- colorRampPalette(c('red4','mistyrose'))

sr_sites <- resp[!is.na(resp$EPT_SppRich),]
sr_sites <- sr_sites[order(sr_sites$EPT_SppRich),]
head(sr_sites)

##convert slope to percent annual change
ave_EPT_SppRich <- mean(metricdata$EPT_SppRich, na.rm=T)
sr_sites$trend_perc <- (sr_sites$EPT_SppRich/ave_EPT_SppRich)*100

#break positive and neg values
sr_neg <- sr_sites$trend_perc[(1:length(sr_sites$trend_perc))[sr_sites$EPT_SppRich <= 0]]
sr_pos <- sr_sites$trend_perc[(1:length(sr_sites$trend_perc))[sr_sites$EPT_SppRich > 0]]

#This adds a column of color values
sr_col_neg <- pal_neg(50)[as.numeric(cut(sr_neg,breaks = 50))]
unique(sr_col_neg)
sr_col_pos <- pal_pos(50)[as.numeric(cut(sr_pos,breaks = 50))]
unique(sr_col_pos)
sr_sites$sr_col <- c(sr_col_neg, sr_col_pos)

newmap <- getMap(resolution = "high")

col <- rep("grey30", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey80")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

plot(newmap,col=col,
bg="paleturquoise4",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(34, 70),
  asp = 1,lwd=wv
)

sr_sites_asc <- sr_sites[order(-sr_sites$trend_perc),]
points(sr_sites$Longitude_X[sr_sites$trend_perc > 0],sr_sites$Latitude_Y[sr_sites$trend_perc > 0],pch = 20,col = alpha(sr_sites$sr_col[sr_sites$trend_perc > 0],0.6),cex=1)

points(sr_sites_asc$Longitude_X[sr_sites_asc$trend_perc <= 0],sr_sites_asc$Latitude_Y[sr_sites_asc$trend_perc <= 0],pch = 20,col = alpha(sr_sites_asc$sr_col[sr_sites_asc$trend_perc <= 0],0.6),cex=1)

lg2 <- c("",round(mean(sr_neg),digits=1),"",0,"",round(mean(sr_pos),digits=1),"")

polygon(x= c(-15,-15,0,0), y= c(60,90,90,60),
        col = "paleturquoise4", border = "paleturquoise4")

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9.15, 72.3,title="",legend=lg2,col =c('midnightblue','white','white','white','white','white','red4'),text.col=c("white"),lty=1,lwd=1,bty="n",cex=0.95)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-12.5,75.2,legend=expression(paste("a) EPT richness (%", y^-1, ")", sep = "")),text.col=c("white"),bty='n',cex=1.7)

############################ insect spp rich #######################################

sr_sites <- resp[!is.na(resp$insect_SppRich),]
sr_sites <- sr_sites[order(sr_sites$insect_SppRich),]
head(sr_sites)

##convert slope to percent annual change
ave_insect_SppRich <- mean(metricdata$insect_SppRich, na.rm=T)
sr_sites$trend_perc <- (sr_sites$insect_SppRich/ave_insect_SppRich)*100

#break positive and neg values
sr_neg <- sr_sites$trend_perc[(1:length(sr_sites$trend_perc))[sr_sites$insect_SppRich <= 0]]
sr_pos <- sr_sites$trend_perc[(1:length(sr_sites$trend_perc))[sr_sites$insect_SppRich > 0]]

#This adds a column of color values
sr_col_neg <- pal_neg(50)[as.numeric(cut(sr_neg,breaks = 50))]
unique(sr_col_neg)
sr_col_pos <- pal_pos(50)[as.numeric(cut(sr_pos,breaks = 50))]
unique(sr_col_pos)
sr_sites$sr_col <- c(sr_col_neg, sr_col_pos)

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

plot(newmap,col=col,
bg="paleturquoise4",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(34, 70),
  asp = 1,lwd=wv
)

sr_sites_asc <- sr_sites[order(-sr_sites$trend_perc),]
points(sr_sites$Longitude_X[sr_sites$trend_perc > 0],sr_sites$Latitude_Y[sr_sites$trend_perc > 0],pch = 20,col = alpha(sr_sites$sr_col[sr_sites$trend_perc > 0],0.6),cex=1)

points(sr_sites_asc$Longitude_X[sr_sites_asc$trend_perc <= 0],sr_sites_asc$Latitude_Y[sr_sites_asc$trend_perc <= 0],pch = 20,col = alpha(sr_sites_asc$sr_col[sr_sites_asc$trend_perc <= 0],0.6),cex=1)

lg2 <- c("",round(mean(sr_neg),digits=1),"",0,"",round(mean(sr_pos),digits=1),"")

polygon(x= c(-15,-15,0,0), y= c(60,90,90,60),
        col = "paleturquoise4", border = "paleturquoise4")

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9.15, 72.3,title="",legend=lg2,col =c('midnightblue','white','white','white','white','white','red4'),text.col=c("white"),lty=1,lwd=1,bty="n",cex=0.95)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-12.5,75.2,legend=expression(paste("b) Insect richness (%", y^-1, ")", sep = "")),text.col=c("white"),bty='n',cex=1.7)

############################ EPT abund #######################################

sr_sites <- resp[!is.na(resp$EPT_Abund),]
sr_sites <- sr_sites[order(sr_sites$EPT_Abund),]
head(sr_sites)

##convert slope to percent annual change
sr_sites$trend_perc <- (10^sr_sites$EPT_Abund-1)*100

#break positive and neg values
sr_neg <- sr_sites$trend_perc[(1:length(sr_sites$trend_perc))[sr_sites$EPT_Abund <= 0]]
sr_pos <- sr_sites$trend_perc[(1:length(sr_sites$trend_perc))[sr_sites$EPT_Abund > 0]]

#This adds a column of color values
sr_col_neg <- pal_neg(50)[as.numeric(cut(sr_neg,breaks = 50))]
unique(sr_col_neg)
sr_col_pos <- pal_pos(50)[as.numeric(cut(sr_pos,breaks = 50))]
unique(sr_col_pos)
sr_sites$sr_col <- c(sr_col_neg, sr_col_pos)

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

plot(newmap,col=col,
bg="paleturquoise4",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(34, 70),
  asp = 1,lwd=wv
)

sr_sites_asc <- sr_sites[order(-sr_sites$trend_perc),]
points(sr_sites$Longitude_X[sr_sites$trend_perc > 0],sr_sites$Latitude_Y[sr_sites$trend_perc > 0],pch = 20,col = alpha(sr_sites$sr_col[sr_sites$trend_perc > 0],0.6),cex=1)

points(sr_sites_asc$Longitude_X[sr_sites_asc$trend_perc <= 0],sr_sites_asc$Latitude_Y[sr_sites_asc$trend_perc <= 0],pch = 20,col = alpha(sr_sites_asc$sr_col[sr_sites_asc$trend_perc <= 0],0.6),cex=1)

lg2 <- c("",round(mean(sr_neg),digits=1),"",0,"",round(mean(sr_pos),digits=1),"")

polygon(x= c(-15,-15,0,0), y= c(60,90,90,60),
        col = "paleturquoise4", border = "paleturquoise4")

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9.15, 72.3,title="",legend=lg2,col =c('midnightblue','white','white','white','white','white','red4'),text.col=c("white"),lty=1,lwd=1,bty="n",cex=0.95)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-12.5,75.2,legend=expression(paste("c) EPT abundance (%", y^-1, ")", sep = "")),text.col=c("white"),bty='n',cex=1.7)

############################ insect abund #######################################

sr_sites <- resp[!is.na(resp$insect_Abund),]
sr_sites <- sr_sites[order(sr_sites$insect_Abund),]
head(sr_sites)

##convert slope to percent annual change
sr_sites$trend_perc <- (10^sr_sites$insect_Abund-1)*100

#break positive and neg values
sr_neg <- sr_sites$trend_perc[(1:length(sr_sites$trend_perc))[sr_sites$insect_Abund <= 0]]
sr_pos <- sr_sites$trend_perc[(1:length(sr_sites$trend_perc))[sr_sites$insect_Abund > 0]]

#This adds a column of color values
sr_col_neg <- pal_neg(50)[as.numeric(cut(sr_neg,breaks = 50))]
unique(sr_col_neg)
sr_col_pos <- pal_pos(50)[as.numeric(cut(sr_pos,breaks = 50))]
unique(sr_col_pos)
sr_sites$sr_col <- c(sr_col_neg, sr_col_pos)

par(mar=c(0,0,0,0))
newmap <- getMap(resolution = "high")

plot(newmap,col=col,
bg="paleturquoise4",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(34, 70),
  asp = 1,lwd=wv
)

sr_sites_asc <- sr_sites[order(-sr_sites$trend_perc),]
points(sr_sites$Longitude_X[sr_sites$trend_perc > 0],sr_sites$Latitude_Y[sr_sites$trend_perc > 0],pch = 20,col = alpha(sr_sites$sr_col[sr_sites$trend_perc > 0],0.6),cex=1)

points(sr_sites_asc$Longitude_X[sr_sites_asc$trend_perc <= 0],sr_sites_asc$Latitude_Y[sr_sites_asc$trend_perc <= 0],pch = 20,col = alpha(sr_sites_asc$sr_col[sr_sites_asc$trend_perc <= 0],0.6),cex=1)

lg2 <- c("",round(mean(sr_neg),digits=1),"",0,"",round(mean(sr_pos),digits=1),"")

polygon(x= c(-15,-15,0,0), y= c(60,90,90,60),
        col = "paleturquoise4", border = "paleturquoise4")

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9.15, 72.3,title="",legend=lg2,col =c('midnightblue','white','white','white','white','white','red4'),text.col=c("white"),lty=1,lwd=1,bty="n",cex=0.95)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=3)
legend(-12.5,75.2,legend=expression(paste("d) Insect abundance (%", y^-1, ")", sep = "")),text.col=c("white"),bty='n',cex=1.7)

dev.off()
#######################
##########################