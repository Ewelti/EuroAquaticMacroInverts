#### libraries #####

library(ggplot2)
library(ggthemes)

#### read in data ####

DATA1_list <- read.csv("raw-data/All_indices_benthicMacroInverts_AllYears.csv", header=T) # change file name according to the time series to be analyzed
DATA2 <- DATA1_list[!is.na(DATA1_list$site_id_wMissing),]


#### do time-series plots ####

#raw data
ggplot(DATA2)+
  geom_point(aes(x=year, y=spp_richness,group=site_id))+
  facet_wrap(~country)+
  theme_few()

#simple linear fits
