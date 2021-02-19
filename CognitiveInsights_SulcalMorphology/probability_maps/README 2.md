## Projected Labels
This pipeline predicts the location of sulci in a given subject based on its coordinates in other subjects.
- MPM_labels.py: Code to build predictive maps in fsaverage
- MPM_run_probability_maps.sh: Run this code to implement the pipeline. It calls MPM_labels.py to build the maps and then projects these maps from the average space back to individual subjects. 
The default threshold is 0.33. This can be manully edited on lines 281 & 398. 

## Requirements
Expects all sulcal labels to be [projected onto fsaverage](https://github.com/cnl-berkeley/lab_scripts/blob/master/freesurfer/label2label.py), and saved in a project directory with the format `${project_dir}/projected_labels/${label_name}/${sub}.?h.${label_name}.label`

## Citations
Bulk of code written by Jacob Miller, with edits from Ben Parker & Willa Voorhies

contact wvoorhies@berkeley.edu with questions

