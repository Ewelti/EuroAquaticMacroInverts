#### libraries #####

library(tidyverse)
library(ggplot2)
library(ggthemes)

#### read in data ####

DATA1_list <- read.csv("raw-data/All_indices_benthicMacroInverts_AllYears.csv", header=T) # change file name according to the time series to be analyzed
DATA2 <- DATA1_list[!is.na(DATA1_list$site_id_wMissing),]

#### cleaning ####

unique(DATA2$month)
DATA2$month[which(DATA2$month == "6to8")] <- "7"
DATA2$month_numeric <- as.numeric(DATA2$month)

#### quick exploratory analysis ####

glimpse(DATA2)

DATA2 %>%
  group_by(country,year) %>%
  summarise(nuYears = length(unique(site_id)))

DATA2 %>%
  group_by(country,site_id) %>%
  summarise(nuYears = length(unique(year)))

DATA2 %>%
  group_by(country,site_id,year,sample_id) %>%
  count()

#### time-series plots ####

#raw data

ggplot(DATA2, aes(x = year, y = spp_richness))+
  geom_point(aes(colour = factor(site_id)))+
  facet_wrap(~ country)+
  theme_few()+
  theme(legend.position = "none")+
  xlab("Year") + ylab("Species richness")

ggsave("plots/raw_timeseries.png", width=10,height=6)

#simple linear fits using geom_smooth

ggplot(DATA2, aes(x=year, y=spp_richness, group=factor(site_id)))+
  facet_wrap(~country)+
  theme_few()+
  geom_smooth(aes(colour = factor(site_id)), method="lm")+
  theme(legend.position = "none")+
  xlab("Year") + ylab("Species richness")

ggsave("plots/simple_linear.png", width=10,height=6)

ggplot(DATA2, aes(x=year, y=spp_richness))+
  facet_wrap(~country)+
  theme_few()+
  geom_smooth(method="lm")+
  xlab("Year") + ylab("Species richness")


#### monthly pattern ?? ####

#amount of data per month
table(DATA2$month)

monthlyData <- DATA2 %>%
  group_by(country,month_numeric) %>%
  summarize(nuYears = length(unique(year)))

ggplot(monthlyData)+
  geom_tile(aes(x = month_numeric, y = country, fill = nuYears))

monthlyData <- DATA2 %>%
  group_by(country,year) %>%
  summarize(month_range = max(month_numeric) - min(month_numeric))
            
ggplot(monthlyData)+
      geom_tile(aes(x = year, y = country, fill = month_range))+
  scale_fill_viridis_c()
            
#seasonal biodiversity patterns
            
ggplot(DATA2, aes(x=month_numeric, y=spp_richness))+
  geom_point(aes(colour = year))+
  facet_wrap(~country)+
  theme_few()+
  scale_colour_viridis_c()+
  xlab("month") + ylab("Species richness")
