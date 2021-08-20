
library(scales)

##Set working directory
##setwd("C:/Users/ewelti/Desktop/aquatic_data")
setwd(dir = "C:/Users/Ellen/Desktop/aquatic_data/")

# attach data
go <- read.csv("All_siteLevel_and_glmOutput.csv", header=T) # change file name according to the time series to be analyzed
attach(go)
head(go)

se <- function(x) sd(x)/sqrt(length(x))

par(mar=c(4,6,0.4,0.4))
##
sp<-SppRich_Est[!is.na(SppRich_Est)]
#average species richness= 27.28712314 species
percChange_perYr<-(sr/27.28712314)*100
d <- density(percChange_perYr_sr)
a <- (max(d$y)+(max(d$y)/10)) *-4
b <- (max(d$y)+(max(d$y)/10))*1
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n')
title(xlab="Percent change", line=2.4,cex.lab=1.3)
axis(2, at=0, labels="Spp R", las=1,cex.axis=1.3)
##
value <- 0
polygon(c(d$x[d$x >= value ], value),
        c(d$y[d$x >= value ], 0),
 	  col = alpha("cornflowerblue", 1),
        border = alpha("cornflowerblue", 1))
polygon(c(d$x[d$x <= value ], value),
        c(d$y[d$x <= value ], 0),
 	  col = alpha("tomato", 1),
        border = alpha("tomato", 1))
me <- mean(percChange_perYr)
ste <- se(percChange_perYr)
sd <- sd(percChange_perYr)
yy <- (4/5*(b-a)+a)
points(x=me, y=yy, lwd=1,pch=20)
arrows(me-ste, yy, me+ste, yy, length=0.05, angle=90, code=3,lwd=3)
arrows(me-sd, yy, me+sd, yy, length=0.05, angle=90, code=3,lwd=3)

##
srr<-SppRichRare_Est[!is.na(SppRichRare_Est)]
#average Rareified species richness= 19.2450259
percChange_perYr<-(srr/19.2450259)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-3
b <- (max(d$y)+(max(d$y)/10))*2
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
axis(2, at=0, labels="Sn", las=1,cex.axis=1.3)
##
value <- 0
polygon(c(d$x[d$x >= value ], value),
        c(d$y[d$x >= value ], 0),
 	  col = alpha("cornflowerblue", 1),
        border = alpha("cornflowerblue", 1))
polygon(c(d$x[d$x <= value ], value),
        c(d$y[d$x <= value ], 0),
 	  col = alpha("tomato", 1),
        border = alpha("tomato", 1))
me <- mean(percChange_perYr)
ste <- se(percChange_perYr)
sd <- sd(percChange_perYr)
yy <- (3/5*(b-a)+a)
points(x=me, y=yy, lwd=1,pch=20)
arrows(me-ste, yy, me+ste, yy, length=0.05, angle=90, code=3,lwd=3)
arrows(me-sd, yy, me+sd, yy, length=0.05, angle=90, code=3,lwd=3)

##
ab<-Abun_Est[!is.na(Abun_Est)]
percChange_perYr <- ab*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-2
b <- (max(d$y)+(max(d$y)/10))*3
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
axis(2, at=0, labels="Abun", las=1,cex.axis=1.3)
##
value <- 0
polygon(c(d$x[d$x >= value ], value),
        c(d$y[d$x >= value ], 0),
 	  col = alpha("cornflowerblue", 1),
        border = alpha("cornflowerblue", 1))
polygon(c(d$x[d$x <= value ], value),
        c(d$y[d$x <= value ], 0),
 	  col = alpha("tomato", 1),
        border = alpha("tomato", 1))
me <- mean(percChange_perYr_sr)
ste <- se(percChange_perYr_sr)
sd <- sd(percChange_perYr_sr)
yy <- (2/5*(b-a)+a)
points(x=me, y=yy, lwd=1,pch=20)
arrows(me-ste, yy, me+ste, yy, length=0.05, angle=90, code=3,lwd=3)
arrows(me-sd, yy, me+sd, yy, length=0.05, angle=90, code=3,lwd=3)

##
spie<-S_PIE_Est[!is.na(S_PIE_Est)]
percChange_perYr <- spie*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10)) *-1
b <- (max(d$y)+(max(d$y)/10))*4
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
axis(2, at=0, labels="S_PIE", las=1,cex.axis=1.3)
##
value <- 0
polygon(c(d$x[d$x >= value ], value),
        c(d$y[d$x >= value ], 0),
 	  col = alpha("cornflowerblue", 1),
        border = alpha("cornflowerblue", 1))
polygon(c(d$x[d$x <= value ], value),
        c(d$y[d$x <= value ], 0),
 	  col = alpha("tomato", 1),
        border = alpha("tomato", 1))
me <- mean(percChange_perYr)
ste <- se(percChange_perYr)
sd <- sd(percChange_perYr)
yy <- (1/5*(b-a)+a)
points(x=me, y=yy, lwd=1,pch=20)
arrows(me-ste, yy, me+ste, yy, length=0.05, angle=90, code=3,lwd=3)
arrows(me-sd, yy, me+sd, yy, length=0.05, angle=90, code=3,lwd=3)

##
to<-TurnO_Est[!is.na(TurnO_Est)]
percChange_perYr <- (to/0.542933401)*100
d <- density(percChange_perYr)
a <- (max(d$y)+(max(d$y)/10))*0
b <- (max(d$y)+(max(d$y)/10))*5
par(new=TRUE)
plot(d, main="",ylab="",xlab="",cex.lab=2,xlim=c(-10,10),ylim=c(a,b),col="white",yaxt='n',xaxt='n')
axis(2, at=0, labels="Turnover", las=1,cex.axis=1.3)
##
value <- 0
polygon(c(d$x[d$x >= value ], value),
        c(d$y[d$x >= value ], 0),
 	  col = alpha("cornflowerblue", 1),
        border = alpha("cornflowerblue", 1))
polygon(c(d$x[d$x <= value ], value),
        c(d$y[d$x <= value ], 0),
 	  col = alpha("tomato", 1),
        border = alpha("tomato", 1))
me <- mean(percChange_perYr)
ste <- se(percChange_perYr)
sd <- sd(percChange_perYr)
yy <- (0/5*(b-a)+a)
points(x=me, y=yy, lwd=1,pch=20)
arrows(me-ste, yy, me+ste, yy, length=0.05, angle=90, code=3,lwd=3)
arrows(me-sd, yy, me+sd, yy, length=0.05, angle=90, code=3,lwd=3)

##
abline(v=0,lwd=2)
##