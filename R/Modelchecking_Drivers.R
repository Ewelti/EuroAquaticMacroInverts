##Set working directory
setwd("C:/Users/ewelti/Desktop/git/EuroAquaticMacroInverts/")

### meta-analysis drivers ####

setwd("outputs/Drivers")
getwd()
library(rstan)
library(brms)
library(loo)

#### spp_richness ####
fit <- readRDS("metaanalysis_drivers_spp_richness.rds")
loo_R2(fit)

#check model
#plot(fit)
sr_loo <- loo(fit, cores = getOption("mc.cores", 1))
sr_loo
sr_parento <- as.list(pareto_k_table(sr_loo))
Count_sr <- rbind(sr_parento[[1]],sr_parento[[2]],sr_parento[[3]],sr_parento[[4]])
colnames(Count_sr) <- "SppRich"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
sr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
sr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
sr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
sr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
sr_fixed <- data.frame(Response="spp_richness", sr_fixed_995[,1:4], sr_fixed_975[,3:4],
                 sr_fixed_95[,3:4],sr_fixed_90[,3:4])
sr_fixed <- cbind(rownames(sr_fixed), data.frame(sr_fixed, row.names=NULL))
colnames(sr_fixed)[1] <- "drivers"

#### spp_rich_rare ####
fit <- readRDS("metaanalysis_drivers_spp_rich_rare.rds")

#check model
#plot(fit)
srr_loo <- loo(fit, cores = getOption("mc.cores", 1))
srr_loo
srr_parento <- as.list(pareto_k_table(srr_loo))
Count_srr <- rbind(srr_parento[[1]],srr_parento[[2]],srr_parento[[3]],srr_parento[[4]])
colnames(Count_srr) <- "SppRichRarefied"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
srr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
srr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
srr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
srr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
srr_fixed <- data.frame(Response="spp_richness_rarefied", srr_fixed_995[,1:4], srr_fixed_975[,3:4],
                  srr_fixed_95[,3:4],srr_fixed_90[,3:4])
srr_fixed <- cbind(rownames(srr_fixed), data.frame(srr_fixed, row.names=NULL))
colnames(srr_fixed)[1] <- "drivers"

#### shannonsH ####
fit <- readRDS("metaanalysis_drivers_shannonsH.rds")

#check model
#plot(fit)
shH_loo <- loo(fit, cores = getOption("mc.cores", 1))
shH_loo
shH_parento <- as.list(pareto_k_table(shH_loo))
Count_shH <- rbind(shH_parento[[1]],shH_parento[[2]],shH_parento[[3]],shH_parento[[4]])
colnames(Count_shH) <- "ShannonsH"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
shH_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
shH_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
shH_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
shH_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
shH_fixed <- data.frame(Response="shannonsH", shH_fixed_995[,1:4], shH_fixed_975[,3:4],
                  shH_fixed_95[,3:4],shH_fixed_90[,3:4])
shH_fixed <- cbind(rownames(shH_fixed), data.frame(shH_fixed, row.names=NULL))
colnames(shH_fixed)[1] <- "drivers"

#### E10 ####
fit <- readRDS("metaanalysis_drivers_E10.rds")

#check model
#plot(fit)
e10_loo <- loo(fit, cores = getOption("mc.cores", 1))
e10_loo
e10_parento <- as.list(pareto_k_table(e10_loo))
Count_e10 <- rbind(e10_parento[[1]],e10_parento[[2]],e10_parento[[3]],e10_parento[[4]])
colnames(Count_e10) <- "E10"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
e10_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
e10_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
e10_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
e10_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
e10_fixed <- data.frame(Response="E10", e10_fixed_995[,1:4], e10_fixed_975[,3:4],
                  e10_fixed_95[,3:4],e10_fixed_90[,3:4])
e10_fixed <- cbind(rownames(e10_fixed), data.frame(e10_fixed, row.names=NULL))
colnames(e10_fixed)[1] <- "drivers"

#### abundance ####
fit <- readRDS("metaanalysis_drivers_abundance.rds")

#check model
#plot(fit)
ab_loo <- loo(fit, cores = getOption("mc.cores", 1))
ab_loo
ab_parento <- as.list(pareto_k_table(ab_loo))
Count_ab <- rbind(ab_parento[[1]],ab_parento[[2]],ab_parento[[3]],ab_parento[[4]])
colnames(Count_ab) <- "Abundance"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
abund_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
abund_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
abund_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
abund_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
abund_fixed <- data.frame(Response="abundance", abund_fixed_995[,1:4], abund_fixed_975[,3:4],
                    abund_fixed_95[,3:4],abund_fixed_90[,3:4])
abund_fixed <- cbind(rownames(abund_fixed), data.frame(abund_fixed, row.names=NULL))
colnames(abund_fixed)[1] <- "drivers"

#### turnover ####
fit <- readRDS("metaanalysis_drivers_turnover.rds")

#check model
#plot(fit)
turn_loo <- loo(fit, cores = getOption("mc.cores", 1))
turn_loo
turn_parento <- as.list(pareto_k_table(turn_loo))
Count_turn <- rbind(turn_parento[[1]],turn_parento[[2]],turn_parento[[3]],turn_parento[[4]])
colnames(Count_turn) <- "turnover"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
turn_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
turn_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
turn_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
turn_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
turn_fixed <- data.frame(Response="turnover", turn_fixed_995[,1:4], turn_fixed_975[,3:4],
                   turn_fixed_95[,3:4],turn_fixed_90[,3:4])
turn_fixed <- cbind(rownames(turn_fixed), data.frame(turn_fixed, row.names=NULL))
colnames(turn_fixed)[1] <- "drivers"

#### F_to ####
fit <- readRDS("metaanalysis_drivers_F_to.rds")

#check model
#plot(fit)
fto_loo <- loo(fit, cores = getOption("mc.cores", 1))
fto_loo
fto_parento <- as.list(pareto_k_table(fto_loo))
Count_fto <- rbind(fto_parento[[1]],fto_parento[[2]],fto_parento[[3]],fto_parento[[4]])
colnames(Count_fto) <- "func_turnover"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
fto_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
fto_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fto_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fto_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fto_fixed <- data.frame(Response="func_turnover", fto_fixed_995[,1:4], fto_fixed_975[,3:4],
                  fto_fixed_95[,3:4],fto_fixed_90[,3:4])
fto_fixed <- cbind(rownames(fto_fixed), data.frame(fto_fixed, row.names=NULL))
colnames(fto_fixed)[1] <- "drivers"

#### FRic ####
fit <- readRDS("metaanalysis_drivers_FRic.rds")

#check model
#plot(fit)
fric_loo <- loo(fit, cores = getOption("mc.cores", 1))
fric_loo
fric_parento <- as.list(pareto_k_table(fric_loo))
Count_fric <- rbind(fric_parento[[1]],fric_parento[[2]],fric_parento[[3]],fric_parento[[4]])
colnames(Count_fric) <- "func_rich"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
fric_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
fric_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fric_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fric_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fric_fixed <- data.frame(Response="func_rich", fric_fixed_995[,1:4], fric_fixed_975[,3:4],
                   fric_fixed_95[,3:4],fric_fixed_90[,3:4])
fric_fixed <- cbind(rownames(fric_fixed), data.frame(fric_fixed, row.names=NULL))
colnames(fric_fixed)[1] <- "drivers"

#### FEve ####
fit <- readRDS("metaanalysis_drivers_FEve.rds")

#check model
#plot(fit)
feve_loo <- loo(fit, cores = getOption("mc.cores", 1))
feve_loo
feve_parento <- as.list(pareto_k_table(feve_loo))
Count_feve <- rbind(feve_parento[[1]],feve_parento[[2]],feve_parento[[3]],feve_parento[[4]])
colnames(Count_feve) <- "func_even"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
feve_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
feve_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
feve_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
feve_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
feve_fixed <- data.frame(Response="func_even", feve_fixed_995[,1:4], feve_fixed_975[,3:4],
                   feve_fixed_95[,3:4],feve_fixed_90[,3:4])
feve_fixed <- cbind(rownames(feve_fixed), data.frame(feve_fixed, row.names=NULL))
colnames(feve_fixed)[1] <- "drivers"

#### FDiv ####
fit <- readRDS("metaanalysis_drivers_FDiv.rds")

#check model
#plot(fit)
fdiv_loo <- loo(fit, cores = getOption("mc.cores", 1))
fdiv_loo
fdiv_parento <- as.list(pareto_k_table(fdiv_loo))
Count_fdiv <- rbind(fdiv_parento[[1]],fdiv_parento[[2]],fdiv_parento[[3]],fdiv_parento[[4]])
colnames(Count_fdiv) <- "func_diverg"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
fdiv_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
fdiv_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
fdiv_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
fdiv_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
fdiv_fixed <- data.frame(Response="func_diverg", fdiv_fixed_995[,1:4], fdiv_fixed_975[,3:4],
                   fdiv_fixed_95[,3:4],fdiv_fixed_90[,3:4])
fdiv_fixed <- cbind(rownames(fdiv_fixed), data.frame(fdiv_fixed, row.names=NULL))
colnames(fdiv_fixed)[1] <- "drivers"

#### RaoQ ####
fit <- readRDS("metaanalysis_drivers_RaoQ.rds")

#check model
#plot(fit)
raoq_loo <- loo(fit, cores = getOption("mc.cores", 1))
raoq_loo
raoq_parento <- as.list(pareto_k_table(raoq_loo))
Count_raoq <- rbind(raoq_parento[[1]],raoq_parento[[2]],raoq_parento[[3]],raoq_parento[[4]])
colnames(Count_raoq) <- "RaoQ"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
raoq_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
raoq_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
raoq_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
raoq_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
raoq_fixed <- data.frame(Response="RaoQ", raoq_fixed_995[,1:4], raoq_fixed_975[,3:4],
                   raoq_fixed_95[,3:4],raoq_fixed_90[,3:4])
raoq_fixed <- cbind(rownames(raoq_fixed), data.frame(raoq_fixed, row.names=NULL))
colnames(raoq_fixed)[1] <- "drivers"

#### alien_SppRich ####
fit <- readRDS("metaanalysis_drivers_alien_SppRich.rds")

#check model
#plot(fit)
aliensr_loo <- loo(fit, cores = getOption("mc.cores", 1))
aliensr_loo
aliensr_parento <- as.list(pareto_k_table(aliensr_loo))
Count_aliensr <- rbind(aliensr_parento[[1]],aliensr_parento[[2]],aliensr_parento[[3]],aliensr_parento[[4]])
colnames(Count_aliensr) <- "alien_sppRich"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
aliensr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
aliensr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
aliensr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
aliensr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
aliensr_fixed <- data.frame(Response="alien_SppRich", aliensr_fixed_995[,1:4], aliensr_fixed_975[,3:4],
                      aliensr_fixed_95[,3:4],aliensr_fixed_90[,3:4])
aliensr_fixed <- cbind(rownames(aliensr_fixed), data.frame(aliensr_fixed, row.names=NULL))
colnames(aliensr_fixed)[1] <- "drivers"

#### alien_Abund ####
fit <- readRDS("metaanalysis_drivers_alien_Abund.rds")

#check model
#plot(fit)
alienab_loo <- loo(fit, cores = getOption("mc.cores", 1))
alienab_loo
alienab_parento <- as.list(pareto_k_table(alienab_loo))
Count_alienab <- rbind(alienab_parento[[1]],alienab_parento[[2]],alienab_parento[[3]],alienab_parento[[4]])
colnames(Count_alienab) <- "alien_abund"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
alienab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
alienab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
alienab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
alienab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
alienab_fixed <- data.frame(Response="alien_abund", alienab_fixed_995[,1:4], alienab_fixed_975[,3:4],
                      alienab_fixed_95[,3:4],alienab_fixed_90[,3:4])
alienab_fixed <- cbind(rownames(alienab_fixed), data.frame(alienab_fixed, row.names=NULL))
colnames(alienab_fixed)[1] <- "drivers"

#### abund_nativeSpp ####
fit <- readRDS("metaanalysis_drivers_abund_nativeSpp.rds")

#check model
#plot(fit)
nativeab_loo <- loo(fit, cores = getOption("mc.cores", 1))
nativeab_parento <- as.list(pareto_k_table(nativeab_loo))
Count_nativeab <- rbind(nativeab_parento[[1]],nativeab_parento[[2]],nativeab_parento[[3]],nativeab_parento[[4]])
colnames(Count_nativeab) <- "native_abund"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
nativeab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
nativeab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
nativeab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
nativeab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
nativeab_fixed <- data.frame(Response="native_abund", nativeab_fixed_995[,1:4], nativeab_fixed_975[,3:4],
                       nativeab_fixed_95[,3:4],nativeab_fixed_90[,3:4])
nativeab_fixed <- cbind(rownames(nativeab_fixed), data.frame(nativeab_fixed, row.names=NULL))
colnames(nativeab_fixed)[1] <- "drivers"

#### SppRich_nativeSpp ####
fit <- readRDS("metaanalysis_drivers_SppRich_nativeSpp.rds")

#check model
#plot(fit)
nativesr_loo <- loo(fit, cores = getOption("mc.cores", 1))
nativesr_parento <- as.list(pareto_k_table(nativesr_loo))
Count_nativesr <- rbind(nativesr_parento[[1]],nativesr_parento[[2]],nativesr_parento[[3]],nativesr_parento[[4]])
colnames(Count_nativesr) <- "native_sppRich"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
nativesr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
nativesr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
nativesr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
nativesr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
nativesr_fixed <- data.frame(Response="native_SppRich", nativesr_fixed_995[,1:4], nativesr_fixed_975[,3:4],
                       nativesr_fixed_95[,3:4],nativesr_fixed_90[,3:4])
nativesr_fixed <- cbind(rownames(nativesr_fixed), data.frame(nativesr_fixed, row.names=NULL))
colnames(nativesr_fixed)[1] <- "drivers"

#### EPT_SppRich ####
fit <- readRDS("metaanalysis_drivers_EPT_SppRich.rds")

#check model
#plot(fit)
EPTsr_loo <- loo(fit, cores = getOption("mc.cores", 1))
EPTsr_parento <- as.list(pareto_k_table(EPTsr_loo))
Count_EPTsr <- rbind(EPTsr_parento[[1]],EPTsr_parento[[2]],EPTsr_parento[[3]],EPTsr_parento[[4]])
colnames(Count_EPTsr) <- "EPT_sppRich"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
EPTsr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
EPTsr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
EPTsr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
EPTsr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
EPTsr_fixed <- data.frame(Response="EPT_SppRich", EPTsr_fixed_995[,1:4], EPTsr_fixed_975[,3:4],
                    EPTsr_fixed_95[,3:4],EPTsr_fixed_90[,3:4])
EPTsr_fixed <- cbind(rownames(EPTsr_fixed), data.frame(EPTsr_fixed, row.names=NULL))
colnames(EPTsr_fixed)[1] <- "drivers"

#### EPT_Abund ####
fit <- readRDS("metaanalysis_drivers_EPT_Abund.rds")

#check model
#plot(fit)
EPTab_loo <- loo(fit, cores = getOption("mc.cores", 1))
EPTab_loo
EPTab_parento <- as.list(pareto_k_table(EPTab_loo))
Count_EPTab <- rbind(EPTab_parento[[1]],EPTab_parento[[2]],EPTab_parento[[3]],EPTab_parento[[4]])
colnames(Count_EPTab) <- "EPT_abund"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
EPTab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
EPTab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
EPTab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
EPTab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
EPTab_fixed <- data.frame(Response="EPT_abund", EPTab_fixed_995[,1:4], EPTab_fixed_975[,3:4],
                    EPTab_fixed_95[,3:4],EPTab_fixed_90[,3:4])
EPTab_fixed <- cbind(rownames(EPTab_fixed), data.frame(EPTab_fixed, row.names=NULL))
colnames(EPTab_fixed)[1] <- "drivers"

#### insect_SppRich ####
fit <- readRDS("metaanalysis_drivers_insect_SppRich.rds")

#check model
#plot(fit)
insectsr_loo <- loo(fit, cores = getOption("mc.cores", 1))
insectsr_parento <- as.list(pareto_k_table(insectsr_loo))
Count_insectsr <- rbind(insectsr_parento[[1]],insectsr_parento[[2]],insectsr_parento[[3]],insectsr_parento[[4]])
colnames(Count_insectsr) <- "insect_sppRich"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
insectsr_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
insectsr_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
insectsr_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
insectsr_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
insectsr_fixed <- data.frame(Response="insect_SppRich", insectsr_fixed_995[,1:4], insectsr_fixed_975[,3:4],
                       insectsr_fixed_95[,3:4],insectsr_fixed_90[,3:4])
insectsr_fixed <- cbind(rownames(insectsr_fixed), data.frame(insectsr_fixed, row.names=NULL))
colnames(insectsr_fixed)[1] <- "drivers"

#### insect_Abund ####
fit <- readRDS("metaanalysis_drivers_insect_Abund.rds")

#check model
#plot(fit)
insectab_loo <- loo(fit, cores = getOption("mc.cores", 1))
insectab_parento <- as.list(pareto_k_table(insectab_loo))
Count_insectab <- rbind(insectab_parento[[1]],insectab_parento[[2]],insectab_parento[[3]],insectab_parento[[4]])
colnames(Count_insectab) <- "insect_abund"
#pp_check(fit, nsamples = 100)

#pull out fixed effects
insectab_fixed_995 <- fixef(fit, probs = c(0.005, 0.995))
insectab_fixed_975 <- fixef(fit, probs = c(0.025, 0.975))
insectab_fixed_95 <- fixef(fit, probs = c(0.05, 0.95))
insectab_fixed_90 <- fixef(fit, probs = c(0.1, 0.9))
insectab_fixed <- data.frame(Response="insect_abund", insectab_fixed_995[,1:4], insectab_fixed_975[,3:4],
                       insectab_fixed_95[,3:4],insectab_fixed_90[,3:4])
insectab_fixed <- cbind(rownames(insectab_fixed), data.frame(insectab_fixed, row.names=NULL))
colnames(insectab_fixed)[1] <- "drivers"

#### assemble all model estimates from Driver meta-analysis models #####

Driver_metaanaly_Ests <- rbind(sr_fixed, srr_fixed, shH_fixed, e10_fixed, abund_fixed, #turn_fixed, 
                           fto_fixed, fric_fixed, feve_fixed, fdiv_fixed, raoq_fixed, aliensr_fixed,
                           alienab_fixed, nativesr_fixed, nativeab_fixed, EPTsr_fixed, EPTab_fixed,
                           insectsr_fixed, insectab_fixed)
write.csv(Driver_metaanaly_Ests, "Driver_metaanaly_Ests.csv")

#### assemble model counts from Parento k diagnostic values from Driver meta-analysis models #####

Driver_metaanaly_parento <- cbind(Count_sr, Count_srr, Count_shH, Count_e10, Count_ab, Count_turn, Count_fto,
                               Count_fric, Count_feve, Count_fdiv, Count_raoq, Count_aliensr,
                              Count_alienab, Count_nativesr, Count_nativeab, Count_EPTsr, Count_EPTab,
                              Count_insectsr, Count_insectab)
rownames(Driver_metaanaly_parento) <- c("good[-Inf, 0.5]","ok[0.5, 0.7]","bad[0.7, 1]","verybad[1, Inf]")
write.csv(Driver_metaanaly_parento, "Driver_meta_parento_ModelCounts.csv")

