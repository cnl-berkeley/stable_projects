#!/bin/bash

selectionfile=/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/fmri_cohort/relmatch_parlabels3/namelist_full.txt

echo $selectionfile
echo ${SGE_TASK_ID}
echo '*'

SUB=$(sed -n "$SGE_TASK_ID"p $selectionfile | awk '{print $1}')

echo $SUB

# currently estimating sa uses relmatch-bp (hard-coded)
# pipeline refers to data to be corrected
pipeline=-bp

cd /home/weiner/shakki/scripts-github/analysis/network_generation/4_regress_spatial_autocorrelation

# calculate 32k fc
python fc_32k-to-32k_wcorr-bp.py $SUB

for hemi in lh rh; do
    # estimate spatial autocorrelation on 32k
    python get_adj_cortex.py $SUB $hemi
    python get_sa-32k.py $SUB $hemi

    # correct manually defined sulcal networks for spatial autocorrelation
    python prune_label_fc.py $SUB $hemi $pipeline nat
    python get_adj_label.py $SUB $hemi nat
    python correct_for_distance.py $SUB $hemi $pipeline nat
done
# combine corrected within- and uncorrected between hemisphere correlations
python hemicombine_label_fc.py $SUB $pipeline nat

### same for probabilistic labels

# correct control analyses using sulcal probability maps
for hemi in lh rh; do
    python prune_label_fc.py $SUB $hemi $pipeline prob
    python get_adj_label.py $SUB $hemi prob
    python correct_for_distance.py $SUB $hemi $pipeline prob
done
python hemicombine_label_fc.py $SUB $pipeline prob

