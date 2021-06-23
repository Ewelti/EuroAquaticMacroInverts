# EuroAquaticMacroInverts
This repository contains code relating to a study of changes over time of European freshwater macroinvertebrates by Haase et al. 

The file "AquaticMacroInvert_metadata" contains information on Data providers, site characteristics, the number of sites from each country, the taxa list for all sites, and the functional trait coverage.

The file "All_indices_benthicMacroInverst_AllYears contains data for each site and year on taxonomic diversity metrics, functional diversity metrics, and environmental drivers which are time series (climate, landuse, and nitrogen inputs).

The file "All_siteLevel_and_glmOutput" contains data for each site on slopes of functional and taxonomic diversity indices and environmental drivers.

The R code "gls_BiodiversityMetrics" was used to calculate slopes of all biodiversity metrics for each site.

Ongoing issues:

--Values of FRic (Functional Richness) and Rao's Q are currently being recaluculated using standardized values within each site.

--Number of dams and dam impact scores are based on all connected dams, but some dams are extremely far from sampling sites. These will be recalucalted using a reduced dataset of dams within a certain area from the sampling location.

--Currently land use values refer to landuse within microbasins the site is located within. Landuse for the larger upstream area is currently being caluculated and will later be added.
