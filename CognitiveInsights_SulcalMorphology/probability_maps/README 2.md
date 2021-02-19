## Probability maps 
Code to generate spatial probability maps. The output thresholded and unthresholded maps for each label are available in *fsaverage_MPM*. 
We provide the code to generate spatial probability maps in fsaverage space and project them into individual subject space. If you would like to use the maps we have generated, you can simply download the fsaverage_MPM directory and project these to an individual subject with the freesurfer command *mri_label2label*. 

**To generate probability maps**
- MPM_labels.py: Code to build spatial probability maps in fsaverage. The default threshold is 0.33. This can be manually edited on lines 281 & 398. 
- MPM_run_probability_maps.sh: Run this code to implement the pipeline. It calls MPM_labels.py to build the maps and then can project these maps from the average space back to individual subjects. 


## Requirements
Expects all sulcal labels to be projected onto fsaverage and saved in a project directory with the format `${project_dir}/projected_labels/${label_name}/${sub}.?h.${label_name}.label`

## Citations
Bulk of code written by Jacob Miller, with edits from Ben Parker & Willa Voorhies

contact wvoorhies@berkeley.edu with questions

