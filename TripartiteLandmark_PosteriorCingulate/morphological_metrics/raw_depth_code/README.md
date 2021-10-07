## Pipeline to extract raw depth metrics from sulcal labels
Computes raw sulcal depth and width for a given labels. 
See [Madan, 2019](https://braininformatics.springeropen.com/articles/10.1186/s40708-019-0098-1) for additional details.

Outputs csvs with the following metrics for all specified label(s), subject(s), and hemisphere(s):
  1. Raw sulcal depth
  2. Raw sulcal width
  
## Requirements
1. Have run the following freesurfer commands on all subjects: [recon-all](https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all) & [recon-all -localGI](https://surfer.nmr.mgh.harvard.edu/fswiki/LGI).
  
2. A subjects list in  a .txt file. 
   - *Note*: you will need to have a separate text file made for each combination of variable sulci for each hemisphere:
     - i.e., those with/without the sps, sspls, and icgs-p in each hemisphere

3. Labels created for a sulcus of interest in the format: subject_dir/subject/label/?h.label_name.label
   - ? = l or r
    
4. Create annotation files for each subject list of variable sulci combinations
   - Use the [mris_label2annot](https://surfer.nmr.mgh.harvard.edu/fswiki/mris_label2annot) freesurfer command
   - Put these annotation files in the /label/ directory for each subject

## USAGE
- Uses *matlab version 2017a*.
- Relevant files that need to be edited:
   - *HCP_wrapper.m* - main file for running code, edit accordingly
   - *calcSulc_load.m* - edit to change the annotation file name accessed
   - *calcSulc.m* - edit to change hemisphere accessed
   - *calcSulc_save.m* - edit so the output file name and column names contain the relevant hemisphere & sulcal labels
   
