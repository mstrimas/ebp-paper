# Best Practices for Using eBird Data

[![Status](https://img.shields.io/badge/Status-in%20prep-red.svg?style=flat-square)]()
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

This repository conducts the analysis and generates the figures for the paper *Best practices for making reliable inferences from citizen science data: case study using eBird to estimate species distributions* ([preprint available on bioRxiv](https://www.biorxiv.org/node/242503.external-links.html)). To reproduce this analysis:

1. Install the latest versions of [R](https://cloud.r-project.org/) and [RStudio](https://www.rstudio.com/products/rstudio/download/#download).
2. Request access to the [eBird Basic Dataset](https://ebird.org/data/download/ebd), then download it and the Sampling Event Data. Uncompress them and place the large text files in a central location.
3. Download this repository by clicking the *Clone or download* button on the top right of this page. Unzip the repository and open the RStudio project.
4. Download the Checklist Calibration Index (CCI) file, after filling out the waiver, and place it in the `data/` subdirectory of the R project. If you do not download this file, the code in the repository will still run, however, CCI will not be included as a predictor in the models.
5. This project use [`packrat`](https://rstudio.github.io/packrat/) to manage package dependencies. When you first open the RStudio project, it should automatically install all the necessary packages required for this analysis. The versions of each package will be consistent with the vesions used in the original analysis. Alternatively, install the dependencies with `remotes::install_github("mstrimas/ebp-paper")`.
6. Run each of the scripts in the top-level directory of the RStudio project in sequence.

## Analysis

The analysis is divided up amongst the following R scripts. The paper uses Wood Thrush as an example species; however, the scripts are set up to use White Ibis and Northern Bobwhite as well. To switch species, change the following in each script: `species <- "Wood Thrush"`.

1. **Setup and data preparation:** these scripts set up the analysis and prepare the data that will be used throughout. They only need to be run once. This steps are optional, all the data resulting from them are in the repository, so the remaining scripts can be run without running these.
  - `00_setup.R`: set up for accessing eBird and MODIS data and prepare a series of GIS layers for the plots. Change the directory for the EBD to wherever you stored the text files. To access MODIS data you'll need to [register for a NASA EarthData account](https://urs.earthdata.nasa.gov/users/new), then put your login and password into this script.
  - `01_ebird-data.R`: extract data from the EBD and prepare it for analysis.
  - `02_habitat-covariates.R`: prepare MODIS landcover based habitat covariates for each checklist.
  - `03_cci-calculation.R`: generate the the CCI file to be used as a predictor. For eBird user privacy, this script is hidden in the repository. Instead, download the prepared CCI data and put it in the `data/` subdirectory.
2. **Best practices:** the following scripts demonstrate the best practices we suggest for modeling encouter rate using random forests, occupancy probability with occupancy models, and relative abundance with GAMs.
  - `04_rf-model.R`: model encounter rate with random forests.
  - `05_occupancy-model.R`: model occupancy probability with occupancy models.
  - `06_gam-count-model.R`: model relative abundance with GAMs.
3. **Bad practices:** the following scripts demostrate the impact that not following the best practices has for each of the three model types.
  - `07_bad_1_rf-model.R`: bad practices for modeling encounter rate with random forests.
  - `07_bad_2_occupancy-model.R`: bad practices for modeling occupancy probability with occupancy models.
  - `07_bad_3_gam-count-model.R`: bad practices for modeling relative abundance with GAMs.
  - `07_bad_4_sample-size-effects.R`: explore the effect of sampling size on model performance, comparing bad and best practice models.
4. **Figures:** code to genate additional figures not part of the main analysis.
 - `08_gobal-effort.R`: generate the global effort map used in the paper.