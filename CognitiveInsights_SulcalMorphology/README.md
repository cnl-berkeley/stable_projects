## Cognitive insights from evolutionarily new brain structures in prefrontal cortex
Analysis pipeline and associated data for the manuscript Voorhies et al. (2020) *Cognitive insights from evolutionarily new brain structures in prefrontal cortex*. The preprint can be accessed [here](https://www.biorxiv.org/content/10.1101/2020.11.07.372805v1).

  For questions or additional data requests please email Willa Voorhies (wvoorhies@berkeley.edu)
  
### Instructions for use ### 
- Access the analysis notebook
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/wvoorhies/CognitiveInsights_SulcalMorphology.git/HEAD)

### Dependencies/Requirements ###
- No downloads are required. Binder will build the enviroment with necessary dependencies to allow the notbook to be run online. 
  - requirements.txt will be used by binder to set up dependencies.
  - install.R installs necessary r-packages (tidyverse)
  - runtime.txt sets R environment. Needed for rpy2
 - *note*: initial environment build can be slow. Depending on your system it can take ~15 mins to setup. This is only true for the first time you launch binder. Binder will output logs during this process, **do not** refresh the page. 
  
### Associated data ###
  - All associated data needed to implement statistical analyses and generate figures is included in the repository.
  - List of data files: 
    - ?_morphology.csv: All subject morphological data (depth, thickness)for each sulcal label in each hemisphere in long format for each sample.
    - ?_Depth.csv: Sulcal depth and behavioral data in wide format for each sample
    - ? thick.csv : cortical thickness and behavioral data in wide format for each sample
    - Plots: model predictions formatted for visualization are saved as .csv files in this dir for later use/plotting in R etc. 
  

