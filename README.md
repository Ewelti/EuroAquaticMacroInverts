## EuroAquaticMacroInverts
This repository contains code relating to the study:

#### **The recovery of European freshwater biodiversity has come to a halt** 

Authors:
Haase, P., D.E. Bowler, N.J. Baker, N. Bonada, S. Domisch, J.R. Garcia Marquez, J. Heino, D. Hering, S.C. Jähnig, A. Schmidt-Kloiber, R. Stubbington, F. Altermatt, M. Álvarez-Cabria, G. Amatulli, D. Angeler, G. Archambaud, I. Arrate Jorrín, T. Aspin, I. Azpiroz, I. Bañares, J. Barquín Ortiz, C.L. Bodin, L. Bonacina, R. Bottarin, M. Cañedo-Argüelles, Z. Csabai, T. Datry, O. Davis, E. de Eyto, A. Dohet, G. Dörflinger, E. Drohan, K.A. Eikland, J. England, T.E. Eriksen, V. Evtimova, M.J. Feio, M. Ferréol, M. Floury, M. Forcellini, M.A. Forio, R. Fornaroli, N. Friberg, J.-F. Fruget, G. Georgieva, P. Goethals, M.A.S. Graça, W. Graf, K.-L. Huttunen, T.C.  Jensen, R.K. Johnson, J.I. Jones, J. Kiesel,, L. Kuglerová, A. Larrañaga, P. Leitner, L. L'Hoste, M.-H. Lizée, A.W. Lorenz, A. Maire, J.A. Manzanos Arnaiz, B. Mckie, A. Millán, D. Monteith, T. Muotka, J.F. Murphy, R. Paavola, P. Paril, F.J. Peñas, F. Pilotto, M. Polasek, J.J. Rasmussen, M. Rubio, D. Sánchez-Fernández, L. Sandin, R.B. Schäfer, A. Scotti, L.Q. Shen, A. Skuja, S. Stoll, M. Straka, H. Timm, V.G. Tyufekchieva, I. Tziortzis, Y. Uzunov, G.H. van der Lee, R. Vannevel, E. Varadinova, G. Várbíró, G. Velle, P.F.M. Verdonschot, R.C.M. Verdonschot, Y. Vidinova, P. Wiberg-Larsen, E.A.R. Welti. 

## Metadata file:

"AquaticMacroInvert_metadata.xls" contains information on data providers, site characteristics, and the number of sites from each country. 

Two tabs within this metadata spreedsheet contain metadata for the two csv files (in the outputs folder):
"All_indices_benthicMacroInverts_AllYears_alienzeros.csv" (metadata tab: metadata_siteyear) and 
"All_siteLevel.csv" (metadata tab: metadata_sites).


## outputs folder:

"All_indices_benthicMacroInverts_AllYears_alienzeros.csv" contains data for each site and year on taxonomic diversity metrics, functional diversity metrics, and environmental drivers which are time series (climate, landuse).

"All_siteLevel.csv" contains data for each site on slopes of functional and taxonomic diversity indices and environmental drivers.

outputs_metaAnaylsis folder - trend meta-analysis model outputs

outputs_movingWindow folder - moving window analysis outputs

outputs_driver folder - driver analysis outputs 

outputs_sensitivity folder - sensitivity analyses outputs including high threshold moving window analysis (HTMW2), one-country removal effects (Jackknife), meta-analysis model comparisions (metaanalysisModelComparison), site level high threshold moving window analysis (siteLevel_HTMW), driver and moving window analyses only for sites with taxonomic resolution to species level(SppLevel_outputs), and taxonomic resolution and seasonality sensitivity analyses outputs (TaxonomicSeason)

TaskIDs folder - task IDs for running models and model estimates

ClimateModel_stan folder - contains outputs of models calculating temperature and precipitation trends


## R folder
### Analyses:

Calculation of trends of biodiversity metrics for each of the 1816 sites

Meta-analysis models to get overall trend for each metric (Trend ~ 1 + (1|study_ID) + (1|Country)

Models with drivers for each metric trend (2 step model)

Moving window analysis

Sensitivity analyses

### HPC analysis: 

Scripts for the analysis based on the HPC to be run as follows:

R/HPC_macroinverts_site_trends.R - calculates trends at the site-level (run on HPC)

R/HPC_combine_site_trends.R - aggregates the site-level data from the previous step (run locally)

R/HPC_Meta-analysis.R - synthesis the site-level data using the output of the previous step (run on HPC)

R/HPC_modelchecking.R - examination of the meta-analysis files of the previous step (run locally)

R/HPC-Meta-analysis_drivers.R - fit driver models to the site-level data (run on HPC)


## plots folder

Online Figures.docx - includes all additional online figures not included in main text or Extended Data

Fig1 folder - contains scripts, data, and icons used to make Figure 1

Fig2_DensityPlots folder - contains scripts and plots related to meta-analysis trend results

Fig3_movingWindow folder - contains scripts and plots related to the moving window analyses

Fig4_drivers folder - contains scripts and plots related to driver analyses

descriptive_plots folder - contains scripts and plots descriping data distributions and correlations

Sensitivity folder - contains scripts and plots for sensitivity analyses regarding one-country removal effects (Jackknife), meta-analysis model comparisions (one stage, two stage weighted and unweighed), effects of sampling years & start year, effects of seasonality, and effects of taxonomic resolution





