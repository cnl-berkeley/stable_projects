## Uncovering a tripartite landmark in posterior cingulate cortex
Analysis pipeline and associated data for the manuscript Willbrand et al. (2022) *Uncovering a tripartite landmark in posterior cingulate cortex*. 
The published manuscript can be accessed at *Science Advances* at this [link](https://www.science.org/doi/10.1126/sciadv.abn9516).

  For questions and/or additional data requests please email Ethan Willbrand (ewillbrand@berkeley.edu) and/or Benjamin Parker (benparker@berkeley.edu)
  
### Instructions for use ### 
  - To run statistical analyses and generate relevant figures open and run the included Rmd file (see *IFRMS_analysis_code.Rmd*).
  - Morphological data was first extracted with FreeSurfer (see **manuscript methods**). A sample function to extract the sulcal metrics of interest from cortical surface reconstructions is included in the */morphological_metrics* directory.
  - Raw depth metrics were also extracted with custom-modified version of a recently developed algorithm building on the FreeSurfer pipeline (see [Madan, 2019](https://braininformatics.springeropen.com/articles/10.1186/s40708-019-0098-1)). Example functions are included under the */morphological_metrics/raw_depth_code* directory.
    - Testing these requires access to a FreeSurfer directory. Please contact Ethan Willbrand (ewillbrand@berkeley.edu) and/or Benjamin Parker (benparker@berkeley.edu) with any questions.
  - Access to and information on the automated sulcal labeling pipeline can be found [here](https://ilwoolyu.github.io/#Software).

### Dependencies/Requirements ###
  - [R](https://www.r-project.org) and [Rstudio](https://www.rstudio.com/products/rstudio/download/) are required to run Rmd file.

  - All necessary R packages are provided at the beginning of the Rmd file. Please download any that are not already downloaded in your library prior to running in order to run the statistical analyses and generate figures.

  
### Associated data ###
  - All associated data needed to implement statistical analyses and generate figures is included in the repository under the */RMD_csvs* directory.
    
  

