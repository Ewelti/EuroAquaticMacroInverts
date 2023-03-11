## EuroAquaticMacroInverts
This repository contains data and code relating to the study:

### **The recovery of European freshwater biodiversity has come to a halt** 

#### Authors:
Haase, P., D.E. Bowler, N.J. Baker, N. Bonada, S. Domisch, J.R. Garcia Marquez, J. Heino, D. Hering, S.C. Jähnig, A. Schmidt-Kloiber, R. Stubbington, F. Altermatt, M. Álvarez-Cabria, G. Amatulli, D. Angeler, G. Archambaud-Suard, I. Arrate Jorrín, T. Aspin, I. Azpiroz, I. Bañares, J. Barquín Ortiz, C.L. Bodin, L. Bonacina, R. Bottarin, M. Cañedo-Argüelles, Z. Csabai, T. Datry, E. de Eyto, A. Dohet, G. Dörflinger, E. Drohan, K.A. Eikland, J. England, T.E. Eriksen, V. Evtimova, M.J. Feio, M. Ferréol, M. Floury, M. Forcellini, M.A. Forio, R. Fornaroli, N. Friberg, J.-F. Fruget, G. Georgieva, P. Goethals, M.A.S. Graça, W. Graf, A. House, K.-L. Huttunen, T.C.  Jensen, R.K. Johnson, J.I. Jones, J. Kiesel,, L. Kuglerová, A. Larrañaga, P. Leitner, L. L'Hoste, M.-H. Lizée, A.W. Lorenz, A. Maire, J.A. Manzanos Arnaiz, B. Mckie, A. Millán, D. Monteith, T. Muotka, J.F. Murphy, D. Ozolins, R. Paavola, P. Paril, F.J. Peñas, F. Pilotto, M. Polasek, J.J. Rasmussen, M. Rubio, D. Sánchez-Fernández, L. Sandin, R.B. Schäfer, A. Scotti, L.Q. Shen, A. Skuja, S. Stoll, M. Straka, H. Timm, V.G. Tyufekchieva, I. Tziortzis, Y. Uzunov, G.H. van der Lee, R. Vannevel, E. Varadinova, G. Várbíró, G. Velle, P.F.M. Verdonschot, R.C.M. Verdonschot, Y. Vidinova, P. Wiberg-Larsen, E.A.R. Welti. 

## Metadata file:

#### _AquaticMacroInvert_metadata.xls_ :
* Contains information on data providers, site characteristics, and the number of sites from each country
* Two tabs within this metadata spreedsheet contain metadata for the two csv files (in the outputs folder):
	* "All_indices_benthicMacroInverts_AllYears_alienzeros.csv" (metadata tab: metadata_siteyear) and 
	* "All_siteLevel.csv" (metadata tab: metadata_sites).


## outputs folder:

#### _All_indices_benthicMacroInverts_AllYears_alienzeros.csv_ : 
Contains data for each site and year on: 
	* taxonomic diversity metrics, 
	* functional diversity metrics, and
	* environmental drivers which are time series (climate, landuse)

#### _All_siteLevel.csv_ :
Contains data for each site on slopes of functional and taxonomic diversity indices and environmental drivers

* _outputs_metaAnaylsis_ folder:
Trend meta-analysis model outputs

* _outputs_movingWindow_ folder:
Moving window analysis outputs

* _outputs_driver_ folder:
Driver analysis outputs 

* _outputs_sensitivity_ folder: 
Sensitivity analyses outputs including: 
	* high threshold moving window analysis (HTMW2), 
	* one-country removal effects (Jackknife), 
	* meta-analysis model comparisions (metaanalysisModelComparison), 
	* site level high threshold moving window analysis (siteLevel_HTMW), 
	* driver and moving window analyses only for sites with taxonomic resolution to species level (SppLevel_outputs), and 
	* taxonomic resolution and seasonality sensitivity analyses outputs (TaxonomicSeason)

* _TaskIDs_ folder:
Task IDs for running models and model estimates

* _ClimateModel_stan_ folder:
Contains outputs of models calculating temperature and precipitation trends


## R folder

* _Initial_Biodiversity_FuncTrait_and_climate_calcs_ folder:
Contains scripts to calculate taxonomic and functional trait metrics for each site year and calculate temperature and precipitation trends

* _stan_models_ folder:
Contains scripts of all stan models used in analyses

* _trend_metaAnalysis_ folder:
Contains the following scripts:
	* HPC_macroinverts_stanmodels: models to calculate trends for each biodiversity metric and each of the 1816 sites
	* HPC_Meta_analysis: Meta-analysis models to synthesize the site-level data using the output of the previous step and get overall trend for each metric (Trend ~ 1 + (1|study_ID) + (1|Country)
	* HPC_modelchecking: examining meta-analysis model fit parameters and calculating probability increasing/decreasing for each biodiversity metric

* _MovingWindow_ folder: 
Contains the following scripts:
	* HPC_MovingWindow_sites: models to calculate trends for each site within each window and biodiversity metric in moving window analysis
	* HPC_Meta_analysis_movingwindow: models to calculate overall estimates of trends within each window and biodiverisity metric
	* MetaMeta_movingwindowEsts_Yr: models to calculate overall linear year effects on biodiversity trajectories in moving window analyses

* _Driver_ folder:
Contains the following scripts:
	* HPC_Meta_analysis_drivers: models to caluculate driver effects on biodiversity metrics using site-level data
	* Modelchecking_Drivers_horseshoePriors: examining driver model fit parameters

* _Sensitivity_ folder:
Contains scripts for sensitivity checking models including: 
	* effects of trends within countries (Country_effects), 
	* high threshold moving window analysis (HTMW), 
	* using a one stage model rather than the two stage meta-analysis model for trend estimates (OneStage_models), 
	* effects of taxonomic resolution for the moving window analysis and claculating the proportion of pos/neg trends/ window (checking_MovingWindow_parameters), 
	* effects of taxonomic resolution and seasonality on trend estimates (HPC_Sensitivity_analysis), and 
	* examing driver model outputs for sites with taxa IDed to species level only (Modelchecking_Drivers_sppLevel_horseshoe)

* _HPC_combine_site_trends_ script:
Concatenates outputs for many analyses including:  
	* biodiverity trends, 
	* meta-analysis, 
	* moving window analysis, 
	* driver analysis, and 
	* sensitivity analyses


## plots folder

#### _Online Figures.docx_ : 
Includes all additional online figures not included in main text or Extended Data

* _Fig1_ folder: 
Contains scripts, data, and icons used to make Figure 1

* _Fig2_DensityPlots_ folder : 
Contains scripts and plots related to meta-analysis trend results

* _Fig3_movingWindow_ folder: 
Contains scripts and plots related to the moving window analyses

* _Fig4_drivers folder_: 
Contains scripts and plots related to driver analyses

* _descriptive_plots_ folder: 
Contains scripts and plots descriping data distributions and correlations

* _Sensitivity_ folder: 
Contains scripts and plots for sensitivity analyses regarding: 
	* one-country removal effects (Jackknife), 
	* meta-analysis model comparisions (one stage, two stage weighted and unweighed), 
	* effects of sampling years & start year, 
	* effects of seasonality, and 
	* effects of taxonomic resolution





