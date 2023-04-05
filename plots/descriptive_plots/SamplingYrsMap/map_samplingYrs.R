##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/")
#

library(RColorBrewer)
library(rworldmap)
library(scales)
library(maps)
x <- map("world", plot=FALSE)
x$names

##import metadata with sites locations
# attach data
sites <- read.csv("outputs/All_siteLevel.csv", header=T) # change file name according to the time series to be analyzed
head(sites)
attach(sites)

################################################################################
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
tiff(filename = "plots/descriptive_plots/SamplingYrsMap/SampleYrs_map.tiff", width = 7, height = 6, units = 'in', res = 600, compression = 'lzw')

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c("gold",'seagreen1','deepskyblue1','purple'))#

sr_sites <- sites[order(sites$Sampling_years),]

sr_sites$sr_col <- rbPal(50)[as.numeric(cut(sr_sites$Sampling_years,breaks = 50))]
unique(sr_sites$sr_col)

par(mar=c(0,0,0,0))

newmap <- getMap(resolution = "high")

col <- rep("grey85", length(newmap@data$NAME))
col[match(ddf$country, newmap@data$NAME)] <- c("black") #or pal instead of c("grey80")

col2 <- rep("#00000000", length(newmap@data$NAME)) #or alpha("white",0) instead of #00000000
col2[match(ddf$country, newmap@data$NAME)] <- c("grey50")

col2[58:59] <- c("black")#change borders for N.Cyprus
wv<-seq(1, 1, length.out=253) 
wv[58:59] <- c(0)

plot(newmap,col=col,
bg="white",border=col2, #bg="lightblue",border="grey70",
  xlim = c(-10, 34),
  ylim = c(35, 70),
  asp = 1,lwd=wv
)

points(sr_sites$Longitude_X,sr_sites$Latitude_Y,pch = 20,col = alpha(sr_sites$sr_col,0.6),cex=1.5)
lg <- round(seq(min(sr_sites$Sampling_years), max(sr_sites$Sampling_years), by=((max(sr_sites$Sampling_years)-min(sr_sites$Sampling_years))/7)),digits=1)
#legend("topleft",title="Species Richness",legend=lg,col =rbPal(8),lty=1,lwd=22,bty="n")

polygon(x= c(-15,-15,0,0), y= c(60,90,90,60),
        col = "white", border = "white")

co_leng<-length(unique(sr_sites$sr_col))
y <- seq(69.8,61.3, by=(-(69.8-61.3)/(co_leng-1)))
x <- rep(-8,co_leng)
x2 <- rep(-6,co_leng)
legend(-9, 72.3,title="",legend=lg,col =rbPal(8),lty=1,lwd=1,bty="n",cex=0.9)
segments(x, y, x2, y, col= unique(sr_sites$sr_col),lwd=6)
legend(-10.2,72.7,legend=c("Sampling years"),bty='n')

dev.off()
#######################
##########################
##########################