# Sulcal depth in prefrontal cortex: A novel predictor of working memory performance
Analysis pipeline and data for the manuscript Yao et al. (2022) Sulcal depth in prefrontal cortex: A novel predictor of working memory performance.

For questions or additional data requests please email Jewelia Yao (jeweliayao@princeton.edu)

### Instructions for use ###
* Access the analysis notebook [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jewelia-yao/PrefrontalSulcalDepth_WorkingMemory/HEAD)


### Dependencies/Requirements ###
* No downloads are required. Binder will build the enviroment with necessary dependencies to allow the notbook to be run online.
  * requirements.txt will be used by binder to set up dependencies.
  * install.R installs necessary r-packages (tidyverse)
  * runtime.txt sets R environment. Needed for rpy2
* Note: initial environment build can be slow. Depending on your system it can take ~15 mins to setup. This is only true for the first time you launch binder. Binder will output logs during this process, do not refresh the page.

### Associated data ###
* All associated data needed to implement statistical analyses and generate figures is included in the repository.

### List of data files: ###
  * sulcal_depth.csv: sulcal depth and behavioral data in wide format for each sample
  * cortical_thickness.csv: cortical thickness and behavioral data in wide format for each sample
  * LHbackwards_scores.csv: behavioral data for the selected model in long format for each subject
  * anova_?.csv: all subject depth data (mean, standard deviation) for each sulcal label in each hemisphere in long format for rm-anovas.
  * demographics: all subject demographic information
  * figures: model predictions formatted for visualization are saved as .csv files in this dir for later use/plotting in R etc.

