# EuroAquaticMacroInverts
This repository contains code relating to a study of changes over time of European freshwater macroinvertebrates by Haase et al. 

Data files:
"AquaticMacroInvert_metadata.xls" contains information on Data providers, site characteristics, the number of sites from each country, the taxa list for all sites, and the functional trait coverage. Two additional tabs provide metadata information for the two csv files "All_indices_benthicMacroInverts_AllYears.csv" (tab: metadata_siteyearsheet) and "All_siteLevel_and_glmOutput.csv" (tab: metadata_sitesheet).
"All_indices_benthicMacroInverts_AllYears.csv" contains data for each site and year on taxonomic diversity metrics, functional diversity metrics, and environmental drivers which are time series (climate, landuse, and nitrogen inputs).
"All_siteLevel_and_glmOutput.csv" contains data for each site on slopes of functional and taxonomic diversity indices and environmental drivers.

Analyses:
-calculate trends of biodiversity metrics for each of the 1816 sites
-meta-analysis models to get overall trend for each metric (Trend ~ 1 + (1|study_ID) + (1|Country)
-models with drivers for each metric trend (2 step model)
-moving window analysis
-models with drivers for each metric each year (1 step model)

HPC analysis: 
Scripts for the analysis based on the HPC to be run as follows
R/HPC_macroinverts_site_trends.R - calculates trends at the site-level (run on HPC)
R/HPC_combine_site_trends.R - aggregates the site-level data from the previous step (run locally)
R/HPC_Meta-analysis.R - synthesis the site-level data using the output of the previous step (run on HPC)
R/HPC_modelchecking.R - examination of the meta-analysis files of the previous step (run locally)
R/HPC-Meta-analysis_drivers.R - fit driver models to the site-level data (run on HPC)


