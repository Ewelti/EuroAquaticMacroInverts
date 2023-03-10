library(tidyverse)
install.packages("janitor")
library(janitor)

#get all samples
modelOutputs <- list.files("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts/outputs/proportions_country", full.names=TRUE) %>%
                        set_names() %>%
                        map_dfr(readRDS, .id="source")

#unpack file names
names <- modelOutputs %>%
                  filter(!duplicated(source)) %>%
                  group_by(source) %>%
                  mutate(Response = strsplit(source, "__")[[1]][2],
                         Country = strsplit(source, "__")[[1]][3],
                         StartYear = strsplit(source, "__")[[1]][4],
                         StartYear = gsub(".RDS","",StartYear)) %>%
                  ungroup() %>%
                  select(source, Response, Country, StartYear)


#add on and summarise
modelSummaries <- modelOutputs %>%
                    inner_join(., names, by="source") %>%
                    janitor::clean_names() %>%
                    group_by(response, start_year, iteration, draw) %>%
                    summarise(propIncrease = mean(b>0)) %>%
                    group_by(response, start_year) %>%
                    summarise(meanProp = mean(propIncrease),
                              lower95Prop = quantile(propIncrease, 0.025),
                              upper95Prop = quantile(propIncrease, 0.975)) %>%
                    ungroup()

modelSummaries %>%
  mutate(start_year = as.numeric(start_year))%>%
ggplot() +
          geom_line(aes(x=start_year, y=meanProp))+
          geom_ribbon(aes(x=start_year, ymin=lower95Prop, ymax=upper95Prop),alpha=0.5)+
          geom_hline(yintercept=0.5, linetype="dashed")+
          facet_wrap(~response)+
          theme_bw()
#############################################################################
ms <- as.data.frame(modelSummaries)
ma <- ms[ which(ms$start_year <2012), ]
ma$MeanYear <- as.numeric(ma$start_year) + 4.5
names(ma)[names(ma) == 'start_year'] <- 'StartYear'

write.csv(ma, "outputs/MovingWindowProportions.csv")
ma <- read.csv("outputs/MovingWindowProportions.csv")

#### import site level estimates
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/EuroAquaticMacroInverts")
ma1 <- readRDS("outputs/outputs_movingWindow/stanTrends_site_level_movingaverages.rds")
ab <- ma1[ma1$Response == "abundance",]
library(data.table)
DT <- data.table(ab)
sitecount_ab <- DT[, .(total_sites = length(unique(site_id))), by = StartYear]

## merge site number
movav1 <- merge(ma,sitecount_ab,by="StartYear", all=F)

#movav1$transnum <- log10(movav1$total_sites)-2.4

##
tiff(filename = "plots/Fig3_movingWindow/proportionPosNeg.tiff", width = 7, height = 6, units = 'in', res = 600, compression = 'lzw')

par(mfrow=c(2,2),mar=c(4,5,0.2,0.2))
count_sign_mw <- movav1[which(movav1$response=='spp_richness'),]
plot(0,0, ylim=c(0.4,0.7),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
box(lwd=2)
title(ylab="Proportion sites w/ + trends", line=3.2,cex.lab=1.1)
title(xlab="Mean year of moving window", line=2.5,cex.lab=1.1)
legend("bottomright", bty="n", legend="a, Taxon richness",cex=1.5)

count_sign_mw <- movav1[which(movav1$response=='abundance'),]
plot(0,0, ylim=c(0.4,0.7),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
box(lwd=2)
title(ylab="Proportion sites w/ + trends", line=3.2,cex.lab=1.1)
title(xlab="Mean year of moving window", line=2.5,cex.lab=1.1)
legend("topright", bty="n", legend="b, Abundance",cex=1.5)

count_sign_mw <- movav1[which(movav1$response=='FRic'),]
plot(0,0, ylim=c(0.4,0.7),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
box(lwd=2)
title(ylab="Proportion sites w/ + trends", line=3.2,cex.lab=1.1)
title(xlab="Mean year of moving window", line=2.5,cex.lab=1.1)
legend("topright", bty="n", legend="c, Functional richness",cex=1.5)

count_sign_mw <- movav1[which(movav1$response=='FRed'),]
plot(0,0, ylim=c(0.4,0.7),xlim= c(1994,2016), ylab="", xlab="",las=1)
polygon(x = c(0.5, 0.5, 2040, 2040), y = c(-100, 0.5, 0.5, -100), col ="coral1", border = NA)
polygon(x = c(0, 0, 2040, 2040), y = c(100, 0.5, 0.5, 100), col ="lightsteelblue1", border = NA)
polygon(x = c(count_sign_mw$MeanYear,rev(count_sign_mw$MeanYear)), y = c(count_sign_mw$lower95Prop,rev(count_sign_mw$upper95Prop)),col="#00000020", border = NA)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, pch=21,bg=1,cex=1.8)
points (count_sign_mw$meanProp ~count_sign_mw$MeanYear, type="l", lwd=2)
box(lwd=2)
title(ylab="Proportion sites w/ + trends", line=3.2,cex.lab=1.1)
title(xlab="Mean year of moving window", line=2.5,cex.lab=1.1)
legend("topright", bty="n", legend="d, Functional redundancy",cex=1.5)

dev.off()
###############################################################################


