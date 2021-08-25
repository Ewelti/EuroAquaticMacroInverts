##Set working directory
##setwd("C:/Users/ewelti/Desktop/aquatic_data")
#setwd(dir = "C:/Users/Ellen/Desktop/aquatic_data/")

# attach data
go <- read.csv("outputs/All_siteLevel_and_glmOutput.csv", header=T) # change file name according to the time series to be analyzed
attach(go)
head(go)

par(mfrow=c(2,3),mar=c(4,5,0.4,0.4))
#density plots
sr<-SppRich_Est[!is.na(SppRich_Est)]
d <- density(sr)
plot(d, main="",ylab="Density",xlab="",cex.lab=2,ylim=c(0,0.85))
#
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Species richness"),cex=1.5,bty="n")
t.test(SppRich_Est)
#average species richness= 27.52824362 species
percChange_perYr_sr<-(mean(sr)/27.52824362)*100
legend("topright",legend=c("+0.92% per year ***"),cex=1.2,bty="n")

par(mar=c(4,2,0.4,0.4))
sh<-ShanH_Est[!is.na(ShanH_Est)]
d <- density(sh)
plot(d, main="",ylab="",xlab="")
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Shannon's H"),cex=1.5,bty="n")
t.test(ShanH_Est)
#average Shannon's H= 1.978251419 
percChange_perYr_sh<-(mean(sh)/1.978251419)*100
legend("topright",legend=c("+0.31% per year ***"),cex=1.2,bty="n")

ev<-EvenJ_Est[!is.na(EvenJ_Est)]
d <- density(ev)
plot(d, main="",ylab="",xlab="")
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Evenness J"),cex=1.5,bty="n")
t.test(EvenJ_Est)
#average Evenness J= 0.619396526
percChange_perYr_ev<-(mean(ev)/0.619396526)*100
legend("topright",legend=c("-0.07% per year, NS"),cex=1.2,bty="n")

par(mar=c(5,5,0.4,0.4))
to<-TurnO_Est[!is.na(TurnO_Est)]
d <- density(to)
plot(d, main="",xlab="",ylab="Density",cex.lab=2,ylim=c(0,45))
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Taxonomic turnover"),cex=1.5,bty="n")
t.test(TurnO_Est)
#average turnover= 0.54294599
percChange_perYr_to<-(mean(to)/0.54294599)*100
legend("topright",legend=c("-0.25% per year ***"),cex=1.2,bty="n")

par(mar=c(5,2,0.4,0.4))
ab<-Adun_Est[!is.na(Adun_Est)]
d <- density(ab)
plot(d, main="",xlab="Site slopes",ylab="",cex.lab=2)
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Abundance"),cex=1.5,bty="n")
t.test(Adun_Est)
#abundance was log transformed so slope is the same value as percent change
mean(ab)*100
legend("topright",legend=c("+0.8% per year ***"),cex=1.2,bty="n")


#######################################
par(mfrow=c(2,3),mar=c(4,5,0.4,0.4))

fr<-FRic_Est[!is.na(FRic_Est)]
d <- density(fr)
plot(d, main="",xlab="",ylab="Density",cex.lab=2,ylim=c(0,18))
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Functional richness"),cex=1.5,bty="n")
t.test(FRic_Est)
#FRic was log transformed so slope is the same value as percent change
mean(fr)*100
legend("topright",legend=c("+1.56% per year ***"),cex=1.2,bty="n")

par(mar=c(4,2,0.4,0.4))
rq<-RaoQ_Est[!is.na(RaoQ_Est)]
d <- density(rq)
plot(d, main="",xlab="",ylab="")
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Rao's Q"),cex=1.5,bty="n")
t.test(RaoQ_Est)
#average Rao's Q= 39.7098828
percChange_perYr_rq<-(mean(rq)/39.7098828)*100
legend("topright",legend=c("+0.07% per year, NS"),cex=1.2,bty="n")

fe<-FEve_Est[!is.na(FEve_Est)]
d <- density(fe)
plot(d, main="",xlab="",ylab="",ylim=c(0,80))
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Functional evenness"),cex=1.5,bty="n")
t.test(FEve_Est)
#average functional evenness= 0.517910888
percChange_perYr_fe<-(mean(fe)/0.517910888)*100
legend("topright",legend=c("-0.24% per year ***"),cex=1.2,bty="n")

par(mar=c(5,5,0.4,0.4))
ft<-F_to_Est[!is.na(F_to_Est)]
d <- density(ft)
plot(d, main="",xlab="",ylab="Density",cex.lab=2,ylim=c(0,130))
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Functional turnover"),cex=1.5,bty="n")
t.test(F_to_Est)
#average functional turnover= 0.086496081
percChange_perYr_ft<-(mean(ft)/0.086496081)*100
legend("topright",legend=c("-0.72% per year ***"),cex=1.2,bty="n")

par(mar=c(5,2,0.4,0.4))
fd<-FDiv_Est[!is.na(FDiv_Est)]
d <- density(fd)
plot(d, main="",xlab="Site slopes",ylab="",cex.lab=2,ylim=c(0,106))
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Functional divergence"),cex=1.5,bty="n")
t.test(FDiv_Est)
#average functional divergence= 0.826257724
percChange_perYr_fd<-(mean(fd)/0.826257724)*100
legend("topright",legend=c("+0.04% per year *"),cex=1.2,bty="n")

##########################################################
######################################################

#Alien
par(mfrow=c(2,2),mar=c(4,5,0.4,0.4))
aa<-AlienAbun_Est[!is.na(AlienAbun_Est)]
d <- density(aa)
plot(d, main="",xlab="Site slopes",ylab="",cex.lab=2,ylim=c(0,35))
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Alien abundance"),cex=1.5,bty="n")
t.test(AlienAbun_Est)
#alien abundance was log transformed so slope is the same value as percent change
mean(aa)*100
legend("topright",legend=c("+0.28% per year, NS"),cex=1.2,bty="n")

asr<-AlienSppRich_Est[!is.na(AlienSppRich_Est)]
d <- density(asr)
plot(d, main="",xlab="Site slopes",ylab="",cex.lab=2,ylim=c(0,30))
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Alien species richness"),cex=1.5,bty="n")
t.test(asr)
#average alien species richness= 1
percChange_perYr_asr<-(mean(asr)/1)*100
legend("topright",legend=c("+2.61% per year ***"),cex=1.2,bty="n")

aap<-AlienAbunProp_Est[!is.na(AlienAbunProp_Est)]
d <- density(aap)
plot(d, main="",xlab="Site slopes",ylab="",cex.lab=2,ylim=c(0,600))
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Proportion alien abundance"),cex=1.5,bty="n")
t.test(aap)
#average proportion alien abundance= 0.088996796
percChange_perYr_aap<-(mean(aap)/0.088996796)*100
legend("topright",legend=c("+2.64% per year ***"),cex=1.2,bty="n")

asrp<-AlienSppRichProp_Est[!is.na(AlienSppRichProp_Est)]
d <- density(asrp)
plot(d, main="",xlab="Site slopes",ylab="",cex.lab=2,ylim=c(0,600))
x_l <- min(which(d$x <= 0))# Define lower limit of colored area
x_ll <- max(which(d$x > 0)) # Define upper limit of colored area
with(d,polygon(x = c(x[c(x_l, x_l:x_ll, x_ll)]),
     y = c(0, y[x_l:x_ll], 0), col = "firebrick3", border="firebrick3"))
x_h <- min(which(d$x >= 0))
x_hh <- max(which(d$x < 100)) 
with(d,polygon(x = c(x[c(x_h, x_h:x_hh, x_hh)]),
     y = c(0, y[x_h:x_hh], 0), col = "dodgerblue", border="dodgerblue"))
legend("topleft",legend=c("Proportion alien species richness"),cex=1.5,bty="n")
t.test(asrp)
#average proportion alien spp rich= 0.049178315
percChange_perYr_asrp<-(mean(asrp)/0.049178315)*100
legend("topright",legend=c("+3.03% per year ***"),cex=1.2,bty="n")
###############################################