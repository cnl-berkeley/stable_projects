#!/bin/bash

selectionfile=/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/fmri_cohort/relmatch_parlabels3/namelist_full.txt

echo $selectionfile
echo ${SGE_TASK_ID}
echo '*'

SUB=$(sed -n "$SGE_TASK_ID"p $selectionfile | awk '{print $1}')
RUNSTR="1,2,3"

echo $SUB
echo $RUNSTR

cd /home/weiner/shakki/scripts-github/analysis/network_generation/3_functional_connectivity_estimation

python fc_seed-to-voxel_wcorr-bp.py $SUB $RUNSTR
python fc_32k-to-32k_wcorr-bp.py $SUB $RUNSTR
#python fc_seed-to-32k_wcorr-bp.py $SUB $RUNSTR

python fc_label-to-label_wcorr.py $SUB $RUNSTR
python fc_label-to-annot_wcorr.py $SUB $RUNSTR
python fc_label-to-label_wcorr-prob.py $SUB $RUNSTR

declare -a LABELS=("ifs" "painfs_any" "pmfs_a" "pmfs_i" "pmfs_p" "sfs_a" "sfs_p" "prts" "lfms" "aalf" "slos1" "sB" "pips" "mTOS" "iTOS" "IPS-PO" "IPS" "cSTS1" "cSTS2" "cSTS3" "aipsJ")
for i in "${LABELS[@]}"
do
    python fc_label-to-annot_wcorr-rotated.py $SUB $RUNSTR "$i" rh
    python fc_label-to-annot_wcorr-rotated.py $SUB $RUNSTR "$i" lh
done

python fc_label-to-label_wcorr-rotated.py $SUB $RUNSTR
 
