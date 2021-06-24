# EuroAquaticMacroInverts
This repository contains code relating to a study of changes over time of European freshwater macroinvertebrates by Haase et al. 

The file "AquaticMacroInvert_metadata.xls" contains information on Data providers, site characteristics, the number of sites from each country, the taxa list for all sites, and the functional trait coverage. Two additional tabs provide metadata information for the two csv files "All_indices_benthicMacroInverst_AllYears.csv" (tab: metadata_siteyearsheet) and "All_siteLevel_and_glmOutput.csv" (tab: metadata_sitesheet).

The file "All_indices_benthicMacroInverst_AllYears.csv" contains data for each site and year on taxonomic diversity metrics, functional diversity metrics, and environmental drivers which are time series (climate, landuse, and nitrogen inputs).

The file "All_siteLevel_and_glmOutput.csv" contains data for each site on slopes of functional and taxonomic diversity indices and environmental drivers.

The R code "gls_BiodiversityMetrics" was used to calculate slopes of all biodiversity metrics for each site.

Ongoing issues:

--Values of FRic (Functional Richness) and Rao's Q are currently being recaluculated using standardized values within each site.

--Currently land-use values refer to land-use within microbasins the site is located within. Land-use for the larger upstream area is currently being caluculated and will later be added.
