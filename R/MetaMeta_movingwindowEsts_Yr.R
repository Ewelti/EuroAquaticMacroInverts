##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")

rs <- readRDS("outputs/stanTrends_site_level.rds")
head(rs)
# attach data
MA <- read.csv("outputs/movingAve_YrEsts.csv")
MA$yr_c <- (MA$StartYear-1989)
head(MA)
MA = subset(MA, select = -c(X) )
MA <- MA[ which(MA$StartYear <2012), ]

##add libraries
#remove.packages(c("StanHeaders", "rstan"))
#install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library("brms")

#define prior
prior1 = c(set_prior("normal(0,3)", class = "Intercept"))

##################abundance
#subset for specific response
resp_est <- subset(MA, Response == "abundance")

fit <- brm(Estimate|se(Est.Error) ~ yr_c, data = resp_est, iter=5000, init = 0,
            chains = 4, prior = prior1,
            control = list(adapt_delta = 0.90,
                           max_treedepth = 12))

#check model
plot(fit)
sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
sr_loo

#pull out fixed effects
sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
sr_fixed <- list(Response="abundance", sr_fixed_995[2,1:4], sr_fixed_975[2,3:4],
                 sr_fixed_95[2,3:4],sr_fixed_90[2,3:4])
ab_fixed <-data.frame(lapply(sr_fixed, function(x) t(data.frame(x))))

##################spp_richness
#subset for specific response
resp_est <- subset(MA, Response == "spp_richness")

fit <- brm(Estimate|se(Est.Error) ~ yr_c, data = resp_est, iter=5000, init = 0,
            chains = 4, prior = prior1,
            control = list(adapt_delta = 0.90,
                           max_treedepth = 12))

#check model
plot(fit)
sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
sr_loo

#pull out fixed effects
sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
sr_fixed <- list(Response="spp_richness", sr_fixed_995[2,1:4], sr_fixed_975[2,3:4],
                 sr_fixed_95[2,3:4],sr_fixed_90[2,3:4])
spprich_fixed <-data.frame(lapply(sr_fixed, function(x) t(data.frame(x))))

##################FRic
#subset for specific response
resp_est <- subset(MA, Response == "FRic")

fit <- brm(Estimate|se(Est.Error) ~ yr_c, data = resp_est, iter=5000, init = 0,
            chains = 4, prior = prior1,
            control = list(adapt_delta = 0.90,
                           max_treedepth = 12))

#check model
plot(fit)
sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
sr_loo

#pull out fixed effects
sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
sr_fixed <- list(Response="FRic", sr_fixed_995[2,1:4], sr_fixed_975[2,3:4],
                 sr_fixed_95[2,3:4],sr_fixed_90[2,3:4])
fric_fixed <-data.frame(lapply(sr_fixed, function(x) t(data.frame(x))))

##################FRed
#subset for specific response
resp_est <- subset(MA, Response == "FRed")

fit <- brm(Estimate|se(Est.Error) ~ yr_c, data = resp_est, iter=5000, init = 0,
            chains = 4, prior = prior1,
            control = list(adapt_delta = 0.90,
                           max_treedepth = 12))

#check model
plot(fit)
sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
sr_loo

#pull out fixed effects
sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
sr_fixed <- list(Response="FRed", sr_fixed_995[2,1:4], sr_fixed_975[2,3:4],
                 sr_fixed_95[2,3:4],sr_fixed_90[2,3:4])
fred_fixed <-data.frame(lapply(sr_fixed, function(x) t(data.frame(x))))

metameta_Ests <- rbind(spprich_fixed, ab_fixed, fric_fixed, fred_fixed)

write.csv(metameta_Ests,"outputs/metameta_MovYrEsts.csv")
#############################################################################
####check percent change in estimates per yr
resp_est <- subset(MA, Response == "spp_richness")
meta_est <- subset(metameta_Ests, Response == "spp_richness")
rm <- mean(resp_est$Estimate)
mm <- meta_est$Estimate
(mm/rm)*100 #### mean percent change in trend est per yr

resp_est <- subset(MA, Response == "abundance")
meta_est <- subset(metameta_Ests, Response == "abundance")
rm <- mean(resp_est$Estimate)
mm <- meta_est$Estimate
(mm/rm)*100 #### mean percent change in trend est per yr

resp_est <- subset(MA, Response == "FRic")
meta_est <- subset(metameta_Ests, Response == "FRic")
rm <- mean(resp_est$Estimate)
mm <- meta_est$Estimate
(mm/rm)*100 #### mean percent change in trend est per yr

resp_est <- subset(MA, Response == "FRed")
meta_est <- subset(metameta_Ests, Response == "FRed")
rm <- mean(resp_est$Estimate)
mm <- meta_est$Estimate
(mm/rm)*100 #### mean percent change in trend est per yr
##############################################################################
metameta_Ests <- read.csv("outputs/metameta_MovYrEsts.csv")
##plot
tiff(filename = "plots/Fig3_movingWindow/metametaEsts.tiff", width = 7, height = 6, units = 'in', res = 600, compression = 'lzw')

par(mar=c(4,1.2,0.4,1), mfrow=c(2,2))

#### spp_richness #####
sub <- subset(metameta_Ests, Response == "spp_richness")
est <- sub$Estimate
yy <- 1
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(sub$Q0.5,0), ylim=c(0.4,1.6),cex=2, bty="n")
segments(x0=0,y0=0,x1=0,y1=4.18,lty=2, lwd=2,col="grey60")
title(xlab="Estimate", line=2.4,cex.lab=1.3)
legend("topleft", legend=("a, Taxon richness"), bty="n", cex=1.3)
mm <- cbind(yy,est)
mm_neg <- as.data.frame(mm[ which(est < 0),])
mm_pos <- as.data.frame(mm[ which(est > 0),])
points(mm_neg[1,] ~ mm_neg[2,], pch="l",cex=12,col="firebrick2")
points(mm_pos[1,] ~ mm_pos[2,], pch="l",cex=12,col="dodgerblue")
yyy1=c(0.975,1.025,1.025,0.975)
yyy2=c(0.95,1.05,1.05,0.95)
yyy3=c(0.9,1.1,1.1,0.9)
polygon(x=c(if (sub$Q2.5<0) {sub$Q2.5} else {0}, if (sub$Q2.5<0) {sub$Q2.5} else {0}, if (sub$Q97.5<0) {sub$Q97.5} else {0}, if (sub$Q97.5<0) {sub$Q97.5} else {0}),y=(yyy1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q5<0) {sub$Q5} else {0}, if (sub$Q5<0) {sub$Q5} else {0}, if (sub$Q95<0) {sub$Q95} else {0}, if (sub$Q95<0) {sub$Q95} else {0}),y=(yyy2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q10<0) {sub$Q10} else {0}, if (sub$Q10<0) {sub$Q10} else {0}, if (sub$Q90<0) {sub$Q90} else {0}, if (sub$Q90<0) {sub$Q90} else {0}),y=(yyy3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q2.5>0) {sub$Q2.5} else {0}, if (sub$Q2.5>0) {sub$Q2.5} else {0}, if (sub$Q97.5>0) {sub$Q97.5} else {0}, if (sub$Q97.5>0) {sub$Q97.5} else {0}),y=(yyy1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$Q5>0) {sub$Q5} else {0}, if (sub$Q5>0) {sub$Q5} else {0}, if (sub$Q95>0) {sub$Q95} else {0}, if (sub$Q95>0) {sub$Q95} else {0}),y=(yyy2), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$Q10>0) {sub$Q10} else {0}, if (sub$Q10>0) {sub$Q10} else {0}, if (sub$Q90>0) {sub$Q90} else {0}, if (sub$Q90>0) {sub$Q90} else {0}),y=(yyy3), col = "dodgerblue", border = NA)
##

#### abundance #####
sub <- subset(metameta_Ests, Response == "abundance")
est <- sub$Estimate
yy <- 1
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(sub$Q0.5,sub$Q99.5), ylim=c(0.4,1.6),cex=2, bty="n")
segments(x0=0,y0=0,x1=0,y1=4.18,lty=2, lwd=2,col="grey60")
title(xlab="Estimate", line=2.4,cex.lab=1.3)
legend("topleft", legend=("b, Abundance"), bty="n", cex=1.3)
mm <- cbind(yy,est)
mm_neg <- as.data.frame(mm[ which(est < 0),])
mm_pos <- as.data.frame(mm[ which(est > 0),])
points(mm_neg[1,] ~ mm_neg[2,], pch="l",cex=12,col="firebrick2")
points(mm_pos[1,] ~ mm_pos[2,], pch="l",cex=12,col="dodgerblue")
yyy1=c(0.975,1.025,1.025,0.975)
yyy2=c(0.95,1.05,1.05,0.95)
yyy3=c(0.9,1.1,1.1,0.9)
polygon(x=c(if (sub$Q2.5<0) {sub$Q2.5} else {0}, if (sub$Q2.5<0) {sub$Q2.5} else {0}, if (sub$Q97.5<0) {sub$Q97.5} else {0}, if (sub$Q97.5<0) {sub$Q97.5} else {0}),y=(yyy1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q5<0) {sub$Q5} else {0}, if (sub$Q5<0) {sub$Q5} else {0}, if (sub$Q95<0) {sub$Q95} else {0}, if (sub$Q95<0) {sub$Q95} else {0}),y=(yyy2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q10<0) {sub$Q10} else {0}, if (sub$Q10<0) {sub$Q10} else {0}, if (sub$Q90<0) {sub$Q90} else {0}, if (sub$Q90<0) {sub$Q90} else {0}),y=(yyy3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q2.5>0) {sub$Q2.5} else {0}, if (sub$Q2.5>0) {sub$Q2.5} else {0}, if (sub$Q97.5>0) {sub$Q97.5} else {0}, if (sub$Q97.5>0) {sub$Q97.5} else {0}),y=(yyy1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$Q5>0) {sub$Q5} else {0}, if (sub$Q5>0) {sub$Q5} else {0}, if (sub$Q95>0) {sub$Q95} else {0}, if (sub$Q95>0) {sub$Q95} else {0}),y=(yyy2), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$Q10>0) {sub$Q10} else {0}, if (sub$Q10>0) {sub$Q10} else {0}, if (sub$Q90>0) {sub$Q90} else {0}, if (sub$Q90>0) {sub$Q90} else {0}),y=(yyy3), col = "dodgerblue", border = NA)
##

#### FRic #####
sub <- subset(metameta_Ests, Response == "FRic")
est <- sub$Estimate
yy <- 1
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(sub$Q0.5,sub$Q99.5), ylim=c(0.4,1.6),cex=2, bty="n")
segments(x0=0,y0=0,x1=0,y1=4.18,lty=2, lwd=2,col="grey60")
title(xlab="Estimate", line=2.4,cex.lab=1.3)
legend("topleft", legend=("c, Functional richness"), bty="n", cex=1.3)
mm <- cbind(yy,est)
mm_neg <- as.data.frame(mm[ which(est < 0),])
mm_pos <- as.data.frame(mm[ which(est > 0),])
points(mm_neg[1,] ~ mm_neg[2,], pch="l",cex=12,col="firebrick2")
points(mm_pos[1,] ~ mm_pos[2,], pch="l",cex=12,col="dodgerblue")
yyy1=c(0.975,1.025,1.025,0.975)
yyy2=c(0.95,1.05,1.05,0.95)
yyy3=c(0.9,1.1,1.1,0.9)
polygon(x=c(if (sub$Q2.5<0) {sub$Q2.5} else {0}, if (sub$Q2.5<0) {sub$Q2.5} else {0}, if (sub$Q97.5<0) {sub$Q97.5} else {0}, if (sub$Q97.5<0) {sub$Q97.5} else {0}),y=(yyy1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q5<0) {sub$Q5} else {0}, if (sub$Q5<0) {sub$Q5} else {0}, if (sub$Q95<0) {sub$Q95} else {0}, if (sub$Q95<0) {sub$Q95} else {0}),y=(yyy2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q10<0) {sub$Q10} else {0}, if (sub$Q10<0) {sub$Q10} else {0}, if (sub$Q90<0) {sub$Q90} else {0}, if (sub$Q90<0) {sub$Q90} else {0}),y=(yyy3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q2.5>0) {sub$Q2.5} else {0}, if (sub$Q2.5>0) {sub$Q2.5} else {0}, if (sub$Q97.5>0) {sub$Q97.5} else {0}, if (sub$Q97.5>0) {sub$Q97.5} else {0}),y=(yyy1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$Q5>0) {sub$Q5} else {0}, if (sub$Q5>0) {sub$Q5} else {0}, if (sub$Q95>0) {sub$Q95} else {0}, if (sub$Q95>0) {sub$Q95} else {0}),y=(yyy2), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$Q10>0) {sub$Q10} else {0}, if (sub$Q10>0) {sub$Q10} else {0}, if (sub$Q90>0) {sub$Q90} else {0}, if (sub$Q90>0) {sub$Q90} else {0}),y=(yyy3), col = "dodgerblue", border = NA)
##

#### FRed #####
sub <- subset(metameta_Ests, Response == "FRed")
est <- sub$Estimate
yy <- 1
plot(yy ~ est, ylab="",xlab="", yaxt="n", las=1, type="n",xlim=c(sub$Q0.5,sub$Q99.5), ylim=c(0.4,1.6),cex=2, bty="n")
segments(x0=0,y0=0,x1=0,y1=4.18,lty=2, lwd=2,col="grey60")
title(xlab="Estimate", line=2.4,cex.lab=1.3)
legend("topleft", legend=("d, Functional redundancy"), bty="n", cex=1.3)
mm <- cbind(yy,est)
mm_neg <- as.data.frame(mm[ which(est < 0),])
mm_pos <- as.data.frame(mm[ which(est > 0),])
points(mm_neg[1,] ~ mm_neg[2,], pch="l",cex=12,col="firebrick2")
points(mm_pos[1,] ~ mm_pos[2,], pch="l",cex=12,col="dodgerblue")
yyy1=c(0.975,1.025,1.025,0.975)
yyy2=c(0.95,1.05,1.05,0.95)
yyy3=c(0.9,1.1,1.1,0.9)
polygon(x=c(if (sub$Q2.5<0) {sub$Q2.5} else {0}, if (sub$Q2.5<0) {sub$Q2.5} else {0}, if (sub$Q97.5<0) {sub$Q97.5} else {0}, if (sub$Q97.5<0) {sub$Q97.5} else {0}),y=(yyy1), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q5<0) {sub$Q5} else {0}, if (sub$Q5<0) {sub$Q5} else {0}, if (sub$Q95<0) {sub$Q95} else {0}, if (sub$Q95<0) {sub$Q95} else {0}),y=(yyy2), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q10<0) {sub$Q10} else {0}, if (sub$Q10<0) {sub$Q10} else {0}, if (sub$Q90<0) {sub$Q90} else {0}, if (sub$Q90<0) {sub$Q90} else {0}),y=(yyy3), col = "firebrick2", border = NA)
polygon(x=c(if (sub$Q2.5>0) {sub$Q2.5} else {0}, if (sub$Q2.5>0) {sub$Q2.5} else {0}, if (sub$Q97.5>0) {sub$Q97.5} else {0}, if (sub$Q97.5>0) {sub$Q97.5} else {0}),y=(yyy1), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$Q5>0) {sub$Q5} else {0}, if (sub$Q5>0) {sub$Q5} else {0}, if (sub$Q95>0) {sub$Q95} else {0}, if (sub$Q95>0) {sub$Q95} else {0}),y=(yyy2), col = "dodgerblue", border = NA)
polygon(x=c(if (sub$Q10>0) {sub$Q10} else {0}, if (sub$Q10>0) {sub$Q10} else {0}, if (sub$Q90>0) {sub$Q90} else {0}, if (sub$Q90>0) {sub$Q90} else {0}),y=(yyy3), col = "dodgerblue", border = NA)
##

##
dev.off()
##


