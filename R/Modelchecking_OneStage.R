#### function to extract posterior distribution of the trends ####

getTrendProbability <- function(fit){
  mySamples <- posterior_samples(fit,pars="b_Intercept")
  data.frame(probIncrease = mean(mySamples>0),probDecrease = mean(mySamples<0))
}

### end of functions ############################################

### one stage models model checking####

#for Ellen
setwd("C:/Users/Ellen/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/onestage")

getwd()
library(rstan)
library(brms)
library(loo)

##################################################################
#### spp_richness ####
fit <- readRDS("onestage_spp_richness.rds")
loo_R2(fit)

#prob of trend
sr_prob <- getTrendProbability(fit)
sr_prob <- data.frame(Response="spp_richness", sr_prob[,1:2])
sr_prob

#check model
plot(fit)
sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
sr_loo
sr_parento <- as.list(pareto_k_table(sr_loo))
Count_sr <- rbind(sr_parento[[1]],sr_parento[[2]],sr_parento[[3]],sr_parento[[4]])
colnames(Count_sr) <- "SppRich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
sr_fixed <- list(Response="spp_richness", sr_fixed_995, sr_fixed_975[3:4],
                 sr_fixed_95[3:4],sr_fixed_90[3:4])
sr_fixed <-data.frame(lapply(sr_fixed, function(x) t(data.frame(x))))
sr_fixed 

#### spp_rich_rare ####
fit <- readRDS("onestage_spp_rich_rare.rds")
summary(fit)
#prob of trend
srr_prob <- getTrendProbability(fit)
srr_prob <- data.frame(Response="spp_richness_rarefied", srr_prob[,1:2])
srr_prob

#check model
plot(fit)
variables(fit)

srr_loo <- loo(fit, cores = getOption("mc.cores", 1))
srr_loo
srr_parento <- as.list(pareto_k_table(srr_loo))
Count_srr <- rbind(srr_parento[[1]],srr_parento[[2]],srr_parento[[3]],srr_parento[[4]])
colnames(Count_srr) <- "SppRichRarefied"
pp_check(fit, nsamples = 100)

#pull out fixed effects
srr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
srr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
srr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
srr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
srr_fixed <- list(Response="spp_richness_rarefied", srr_fixed_995, srr_fixed_975[3:4],
                  srr_fixed_95[3:4],srr_fixed_90[3:4])
srr_fixed <-data.frame(lapply(srr_fixed, function(x) t(data.frame(x))))
srr_fixed

#### shannonsH ####
fit <- readRDS("onestage_shannonsH.rds")

#prob of trend
shH_prob <- getTrendProbability(fit)
shH_prob <- data.frame(Response="shannonsH", shH_prob[,1:2])
shH_prob

#check model
plot(fit)
shH_loo <- loo(fit, cores = getOption("mc.cores", 1))
save_pars(fit)
shH_loo
shH_parento <- as.list(pareto_k_table(shH_loo))

Count_shH <- rbind(shH_parento[[1]],shH_parento[[2]],shH_parento[[3]],shH_parento[[4]])
colnames(Count_shH) <- "ShannonsH"
pp_check(fit, nsamples = 100)

#pull out fixed effects
shH_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
shH_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
shH_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
shH_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
shH_fixed <- list(Response="shannonsH", shH_fixed_995, shH_fixed_975[3:4],
                  shH_fixed_95[3:4],shH_fixed_90[3:4])
shH_fixed <-data.frame(lapply(shH_fixed, function(x) t(data.frame(x))))
shH_fixed

#### E10 ####
fit <- readRDS("onestage_E10.rds")

#prob of trend
e10_prob <- getTrendProbability(fit)
e10_prob <- data.frame(Response="E10", e10_prob[,1:2])
e10_prob

#check model
plot(fit)
e10_loo <- loo(fit, cores = getOption("mc.cores", 1))
e10_loo
e10_parento <- as.list(pareto_k_table(e10_loo))
Count_e10 <- rbind(e10_parento[[1]],e10_parento[[2]],e10_parento[[3]],e10_parento[[4]])
colnames(Count_e10) <- "E10"
pp_check(fit, nsamples = 100)

#pull out fixed effects
e10_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
e10_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
e10_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
e10_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
e10_fixed <- list(Response="E10", e10_fixed_995, e10_fixed_975[3:4],
                  e10_fixed_95[3:4],e10_fixed_90[3:4])
e10_fixed <-data.frame(lapply(e10_fixed, function(x) t(data.frame(x))))
e10_fixed

#### abundance ####
fit <- readRDS("onestage_abundance.rds")

#prob of trend
ab_prob <- getTrendProbability(fit)
ab_prob <- data.frame(Response="abundance", ab_prob[,1:2])
ab_prob

#check model
plot(fit)
ab_loo <- loo(fit, cores = getOption("mc.cores", 1))
ab_loo
ab_parento <- as.list(pareto_k_table(ab_loo))
Count_ab <- rbind(ab_parento[[1]],ab_parento[[2]],ab_parento[[3]],ab_parento[[4]])
colnames(Count_ab) <- "Abundance"
pp_check(fit, nsamples = 100)

#pull out fixed effects
abund_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
abund_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
abund_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
abund_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
abund_fixed <- list(Response="abundance", abund_fixed_995, abund_fixed_975[3:4],
                    abund_fixed_95[3:4],abund_fixed_90[3:4])
abund_fixed <-data.frame(lapply(abund_fixed, function(x) t(data.frame(x))))
abund_fixed

#### turnover ####
fit <- readRDS("onestage_turnover.rds")

#prob of trend
turn_prob <- getTrendProbability(fit)
turn_prob <- data.frame(Response="turnover", turn_prob[,1:2])
turn_prob

#check model
plot(fit)
turn_loo <- loo(fit, cores = getOption("mc.cores", 1))
turn_loo
turn_parento <- as.list(pareto_k_table(turn_loo))
Count_turn <- rbind(turn_parento[[1]],turn_parento[[2]],turn_parento[[3]],turn_parento[[4]])
colnames(Count_turn) <- "turnover"
pp_check(fit, nsamples = 100)

#pull out fixed effects
turn_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
turn_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
turn_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
turn_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
turn_fixed <- list(Response="turnover", turn_fixed_995, turn_fixed_975[3:4],
                   turn_fixed_95[3:4],turn_fixed_90[3:4])
turn_fixed <-data.frame(lapply(turn_fixed, function(x) t(data.frame(x))))
turn_fixed

#### F_to ####
fit <- readRDS("onestage_F_to.rds")

#prob of trend
fto_prob <- getTrendProbability(fit)
fto_prob <- data.frame(Response="func_turnover", fto_prob[,1:2])
fto_prob

#check model
plot(fit)# looks okay even unweighted now
fto_loo <- loo(fit, cores = getOption("mc.cores", 1))
fto_loo
fto_parento <- as.list(pareto_k_table(fto_loo))
Count_fto <- rbind(fto_parento[[1]],fto_parento[[2]],fto_parento[[3]],fto_parento[[4]])
colnames(Count_fto) <- "func_turnover"
pp_check(fit, nsamples = 100)

#pull out fixed effects
fto_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
fto_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
fto_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
fto_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
fto_fixed <- list(Response="func_turnover", fto_fixed_995, fto_fixed_975[3:4],
                  fto_fixed_95[3:4],fto_fixed_90[3:4])
fto_fixed <-data.frame(lapply(fto_fixed, function(x) t(data.frame(x))))
fto_fixed

#### FRic ####
fit <- readRDS("onestage_FRic.rds")

#prob of trend
fric_prob <- getTrendProbability(fit)
fric_prob <- data.frame(Response="func_rich", fric_prob[,1:2])
fric_prob

#check model
plot(fit)
fric_loo <- loo(fit, cores = getOption("mc.cores", 1))
fric_loo
fric_parento <- as.list(pareto_k_table(fric_loo))
Count_fric <- rbind(fric_parento[[1]],fric_parento[[2]],fric_parento[[3]],fric_parento[[4]])
colnames(Count_fric) <- "func_rich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
fric_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
fric_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
fric_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
fric_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
fric_fixed <- list(Response="func_rich", fric_fixed_995, fric_fixed_975[3:4],
                   fric_fixed_95[3:4],fric_fixed_90[3:4])
fric_fixed <-data.frame(lapply(fric_fixed, function(x) t(data.frame(x))))
fric_fixed

#### FEve ####
fit <- readRDS("onestage_FEve.rds")

#prob of trend
feve_prob <- getTrendProbability(fit)
feve_prob <- data.frame(Response="func_even", feve_prob[,1:2])
feve_prob

#check model
plot(fit)
feve_loo <- loo(fit, cores = getOption("mc.cores", 1))
feve_loo
feve_parento <- as.list(pareto_k_table(feve_loo))
Count_feve <- rbind(feve_parento[[1]],feve_parento[[2]],feve_parento[[3]],feve_parento[[4]])
colnames(Count_feve) <- "func_even"
pp_check(fit, nsamples = 100)

#pull out fixed effects
feve_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
feve_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
feve_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
feve_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
feve_fixed <- list(Response="func_even", feve_fixed_995, feve_fixed_975[3:4],
                   feve_fixed_95[3:4],feve_fixed_90[3:4])
feve_fixed <-data.frame(lapply(feve_fixed, function(x) t(data.frame(x))))
feve_fixed

#### FDiv ####
fit <- readRDS("onestage_FDiv.rds")

#prob of trend
fdiv_prob <- getTrendProbability(fit)
fdiv_prob <- data.frame(Response="func_diverg", fdiv_prob[,1:2])
fdiv_prob

#check model
plot(fit)
fdiv_loo <- loo(fit, cores = getOption("mc.cores", 1))
fdiv_loo
fdiv_parento <- as.list(pareto_k_table(fdiv_loo))
Count_fdiv <- rbind(fdiv_parento[[1]],fdiv_parento[[2]],fdiv_parento[[3]],fdiv_parento[[4]])
colnames(Count_fdiv) <- "func_diverg"
pp_check(fit, nsamples = 100)

#pull out fixed effects
fdiv_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
fdiv_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
fdiv_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
fdiv_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
fdiv_fixed <- list(Response="func_diverg", fdiv_fixed_995, fdiv_fixed_975[3:4],
                   fdiv_fixed_95[3:4],fdiv_fixed_90[3:4])
fdiv_fixed <-data.frame(lapply(fdiv_fixed, function(x) t(data.frame(x))))
fdiv_fixed

#### RaoQ ####
fit <- readRDS("onestage_RaoQ.rds")

#prob of trend
raoq_prob <- getTrendProbability(fit)
raoq_prob <- data.frame(Response="RaoQ", raoq_prob[,1:2])
raoq_prob

#check model
plot(fit)
raoq_loo <- loo(fit, cores = getOption("mc.cores", 1))
raoq_loo
raoq_parento <- as.list(pareto_k_table(raoq_loo))
Count_raoq <- rbind(raoq_parento[[1]],raoq_parento[[2]],raoq_parento[[3]],raoq_parento[[4]])
colnames(Count_raoq) <- "RaoQ"
pp_check(fit, nsamples = 100)

#pull out fixed effects
raoq_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
raoq_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
raoq_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
raoq_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
raoq_fixed <- list(Response="RaoQ", raoq_fixed_995, raoq_fixed_975[3:4],
                   raoq_fixed_95[3:4],raoq_fixed_90[3:4])
raoq_fixed <-data.frame(lapply(raoq_fixed, function(x) t(data.frame(x))))
raoq_fixed

#### FRed ####
fit <- readRDS("onestage_FRed.rds")

#prob of trend
FRed_prob <- getTrendProbability(fit)
FRed_prob <- data.frame(Response="FRed", FRed_prob[,1:2])
FRed_prob

#check model
plot(fit)
FRed_loo <- loo(fit, cores = getOption("mc.cores", 1))
FRed_loo
FRed_parento <- as.list(pareto_k_table(FRed_loo))
Count_FRed <- rbind(FRed_parento[[1]],FRed_parento[[2]],FRed_parento[[3]],FRed_parento[[4]])
colnames(Count_FRed) <- "FRed"
pp_check(fit, nsamples = 100)

#pull out fixed effects
FRed_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
FRed_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
FRed_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
FRed_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
FRed_fixed <- list(Response="FRed", FRed_fixed_995, FRed_fixed_975[3:4],
                   FRed_fixed_95[3:4],FRed_fixed_90[3:4])
FRed_fixed <-data.frame(lapply(FRed_fixed, function(x) t(data.frame(x))))
FRed_fixed

#### alien_SppRich ####
fit <- readRDS("onestage_alien_SppRich.rds")

#prob of trend
aliensr_prob <- getTrendProbability(fit)
aliensr_prob <- data.frame(Response="alien_sppRich", aliensr_prob[,1:2])
aliensr_prob

#check model
plot(fit)#bad!! - better unweighted
aliensr_loo <- loo(fit, cores = getOption("mc.cores", 1))
aliensr_loo
aliensr_parento <- as.list(pareto_k_table(aliensr_loo))
Count_aliensr <- rbind(aliensr_parento[[1]],aliensr_parento[[2]],aliensr_parento[[3]],aliensr_parento[[4]])
colnames(Count_aliensr) <- "alien_sppRich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
aliensr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
aliensr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
aliensr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
aliensr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
aliensr_fixed <- list(Response="alien_SppRich", aliensr_fixed_995, aliensr_fixed_975[3:4],
                      aliensr_fixed_95[3:4],aliensr_fixed_90[3:4])
aliensr_fixed <-data.frame(lapply(aliensr_fixed, function(x) t(data.frame(x))))
aliensr_fixed

#### alien_Abund ####
fit <- readRDS("onestage_alien_Abund.rds")

#prob of trend
alienab_prob <- getTrendProbability(fit)
alienab_prob <- data.frame(Response="alien_abund", alienab_prob[,1:2])
alienab_prob

#check model
plot(fit) 
alienab_loo <- loo(fit, cores = getOption("mc.cores", 1))
alienab_loo
alienab_parento <- as.list(pareto_k_table(alienab_loo))
Count_alienab <- rbind(alienab_parento[[1]],alienab_parento[[2]],alienab_parento[[3]],alienab_parento[[4]])
colnames(Count_alienab) <- "alien_abund"
pp_check(fit, nsamples = 100)

#pull out fixed effects
alienab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
alienab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
alienab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
alienab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
alienab_fixed <- list(Response="alien_abund", alienab_fixed_995, alienab_fixed_975[3:4],
                      alienab_fixed_95[3:4],alienab_fixed_90[3:4])
alienab_fixed <-data.frame(lapply(alienab_fixed, function(x) t(data.frame(x))))
alienab_fixed

#### abund_nativeSpp ####
fit <- readRDS("onestage_abund_nativeSpp.rds")

#prob of trend
nativeab_prob <- getTrendProbability(fit)
nativeab_prob <- data.frame(Response="native_abund", nativeab_prob[,1:2])
nativeab_prob

#check model
plot(fit)# seems okay weighted now
nativeab_loo <- loo(fit, cores = getOption("mc.cores", 1))
nativeab_loo
nativeab_parento <- as.list(pareto_k_table(nativeab_loo))
Count_nativeab <- rbind(nativeab_parento[[1]],nativeab_parento[[2]],nativeab_parento[[3]],nativeab_parento[[4]])
colnames(Count_nativeab) <- "native_abund"
pp_check(fit, nsamples = 100)

#pull out fixed effects
nativeab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
nativeab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
nativeab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
nativeab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
nativeab_fixed <- list(Response="native_abund", nativeab_fixed_995, nativeab_fixed_975[3:4],
                       nativeab_fixed_95[3:4],nativeab_fixed_90[3:4])
nativeab_fixed <-data.frame(lapply(nativeab_fixed, function(x) t(data.frame(x))))
nativeab_fixed

#### SppRich_nativeSpp ####
fit <- readRDS("onestage_SppRich_nativeSpp.rds")

#prob of trend
nativesr_prob <- getTrendProbability(fit)
nativesr_prob <- data.frame(Response="native_sppRich", nativesr_prob[,1:2])
nativesr_prob

#check model
plot(fit)
nativesr_loo <- loo(fit, cores = getOption("mc.cores", 1))
nativesr_loo
nativesr_parento <- as.list(pareto_k_table(nativesr_loo))
Count_nativesr <- rbind(nativesr_parento[[1]],nativesr_parento[[2]],nativesr_parento[[3]],nativesr_parento[[4]])
colnames(Count_nativesr) <- "native_sppRich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
nativesr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
nativesr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
nativesr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
nativesr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
nativesr_fixed <- list(Response="native_SppRich", nativesr_fixed_995, nativesr_fixed_975[3:4],
                       nativesr_fixed_95[3:4],nativesr_fixed_90[3:4])
nativesr_fixed <-data.frame(lapply(nativesr_fixed, function(x) t(data.frame(x))))
nativesr_fixed

#### EPT_SppRich ####
fit <- readRDS("onestage_EPT_SppRich.rds")

#prob of trend
EPTsr_prob <- getTrendProbability(fit)
EPTsr_prob <- data.frame(Response="EPT_sppRich", EPTsr_prob[,1:2])
EPTsr_prob

#check model
plot(fit)#bad!! better unweighted
EPTsr_loo <- loo(fit, cores = getOption("mc.cores", 1))
EPTsr_loo
EPTsr_parento <- as.list(pareto_k_table(EPTsr_loo))
Count_EPTsr <- rbind(EPTsr_parento[[1]],EPTsr_parento[[2]],EPTsr_parento[[3]],EPTsr_parento[[4]])
colnames(Count_EPTsr) <- "EPT_sppRich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
EPTsr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
EPTsr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
EPTsr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
EPTsr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
EPTsr_fixed <- list(Response="EPT_SppRich", EPTsr_fixed_995, EPTsr_fixed_975[3:4],
                    EPTsr_fixed_95[3:4],EPTsr_fixed_90[3:4])
EPTsr_fixed <-data.frame(lapply(EPTsr_fixed, function(x) t(data.frame(x))))
EPTsr_fixed

#### EPT_Abund ####
fit <- readRDS("onestage_EPT_Abund.rds")

#prob of trend
EPTab_prob <- getTrendProbability(fit)
EPTab_prob <- data.frame(Response="EPT_abund", EPTab_prob[,1:2])
EPTab_prob

#check model
plot(fit)#bad!! better unweighted
EPTab_loo <- loo(fit, cores = getOption("mc.cores", 1))
EPTab_loo
EPTab_parento <- as.list(pareto_k_table(EPTab_loo))
Count_EPTab <- rbind(EPTab_parento[[1]],EPTab_parento[[2]],EPTab_parento[[3]],EPTab_parento[[4]])
colnames(Count_EPTab) <- "EPT_abund"
pp_check(fit, nsamples = 100)

#pull out fixed effects
EPTab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
EPTab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
EPTab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
EPTab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
EPTab_fixed <- list(Response="EPT_abund", EPTab_fixed_995, EPTab_fixed_975[3:4],
                    EPTab_fixed_95[3:4],EPTab_fixed_90[3:4])
EPTab_fixed <-data.frame(lapply(EPTab_fixed, function(x) t(data.frame(x))))
EPTab_fixed

#### insect_SppRich ####
fit <- readRDS("onestage_insect_SppRich.rds")

#prob of trend
insectsr_prob <- getTrendProbability(fit)
insectsr_prob <- data.frame(Response="insect_sppRich", insectsr_prob[,1:2])
insectsr_prob

#check model
plot(fit)
insectsr_loo <- loo(fit, cores = getOption("mc.cores", 1))
insectsr_loo
insectsr_parento <- as.list(pareto_k_table(insectsr_loo))
Count_insectsr <- rbind(insectsr_parento[[1]],insectsr_parento[[2]],insectsr_parento[[3]],insectsr_parento[[4]])
colnames(Count_insectsr) <- "insect_sppRich"
pp_check(fit, nsamples = 100)

#pull out fixed effects
insectsr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
insectsr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
insectsr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
insectsr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
insectsr_fixed <- list(Response="insect_SppRich", insectsr_fixed_995, insectsr_fixed_975[3:4],
                       insectsr_fixed_95[3:4],insectsr_fixed_90[3:4])
insectsr_fixed <-data.frame(lapply(insectsr_fixed, function(x) t(data.frame(x))))
insectsr_fixed

#### insect_Abund ####
fit <- readRDS("onestage_insect_Abund.rds")

#prob of trend
insectab_prob <- getTrendProbability(fit)
insectab_prob <- data.frame(Response="insect_abund", insectab_prob[,1:2])
insectab_prob

#check model
plot(fit)
insectab_loo <- loo(fit, cores = getOption("mc.cores", 1))
insectab_loo
insectab_parento <- as.list(pareto_k_table(insectab_loo))
Count_insectab <- rbind(insectab_parento[[1]],insectab_parento[[2]],insectab_parento[[3]],insectab_parento[[4]])
colnames(Count_insectab) <- "insect_abund"
pp_check(fit, nsamples = 100)

#pull out fixed effects
insectab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))[2,]
insectab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))[2,]
insectab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))[2,]
insectab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))[2,]
insectab_fixed <- list(Response="insect_abund", insectab_fixed_995, insectab_fixed_975[3:4],
                       insectab_fixed_95[3:4],insectab_fixed_90[3:4])
insectab_fixed <-data.frame(lapply(insectab_fixed, function(x) t(data.frame(x))))
insectab_fixed

#### assemble all model estimates from meta-analysis models #####

Yr_OneStage_Ests <- rbind(sr_fixed, srr_fixed, shH_fixed, e10_fixed, abund_fixed, turn_fixed, 
                           fto_fixed, fric_fixed, feve_fixed, fdiv_fixed, raoq_fixed, FRed_fixed, aliensr_fixed,
                           alienab_fixed, nativesr_fixed, nativeab_fixed, EPTsr_fixed, EPTab_fixed,
                           insectsr_fixed, insectab_fixed)
write.csv(Yr_OneStage_Ests, "Yr_OneStage_Ests.csv")

#### assemble all probabilities of increases/decreases from meta-analysis models #####

Yr_OneStage_probs <- rbind(sr_prob, srr_prob, shH_prob, e10_prob, ab_prob, turn_prob, 
                            fto_prob, fric_prob, feve_prob, fdiv_prob, raoq_prob, FRed_prob, aliensr_prob,
                            alienab_prob, nativesr_prob, nativeab_prob, EPTsr_prob, EPTab_prob,
                            insectsr_prob, insectab_prob)
write.csv(Yr_OneStage_probs, "Yr_OneStage_probabilities.csv")

#### assemble model counts from Parento k diagnostic values from meta-analysis models #####

Yr_OneStage_parento <- cbind(Count_sr, Count_srr, Count_shH, Count_e10, Count_ab, Count_turn, Count_fto,
                              Count_fric, Count_feve, Count_fdiv, Count_raoq, Count_FRed, Count_aliensr,
                              Count_alienab, Count_nativesr, Count_nativeab, Count_EPTsr, Count_EPTab,
                              Count_insectsr, Count_insectab)
rownames(Yr_OneStage_parento) <- c("good[-Inf, 0.5]","ok[0.5, 0.7]","bad[0.7, 1]","verybad[1, Inf]")
write.csv(Yr_OneStage_parento, "Yr_OneStage_parento_ModelCounts.csv")