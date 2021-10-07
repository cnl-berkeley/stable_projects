## Pipeline to extract morphological metrics from sulcal labels
Pulls common morphological metrics for given Freesurfer labels.
Computes mean/max/normalized sulcal depth for a given labels. 
Normalized mean sulcal depth for a label is given as a percentage of maximum depth in a hemisphere. 

Outputs a single CSV with the following metrics for all given labels and subjects:
  1. Mean/sd Cortical Thickness
  2. Surface Area
  3. Gray matter volume 
  4. Mean/Gaussian/intrinsic curvature
  5.Folding index
  6. mean/max/normalized sulcal depth
  
## Requirements
1. have run the following freesurfer commands on all subjects: recon-all, mris_anatomicals: output in nested in a "label_stats" directory in subjects label folder.
  (ie. ../<label>/<label_stats>/<hemi>.<label>.label.stats.txt)
  
2. a subjects list in  a .txt file.

3. Labels created for a sulcus of interest in the format: subject_dir/subject/<label>/
    
4. nibabel and nilearn must be installed and operational

## USAGE
sulcal_morph.py '<subject_dir>' '<sublist.txt>' <labels> '<out_dir>'
- subjects_dir = path to freesurfer subjects directory 
- sublist = path to subjects list text file
- labels = labels of interest as a list of stirngs eg. ['POS', 'MCGS', 'sbps']
- outdir = path to output directory
