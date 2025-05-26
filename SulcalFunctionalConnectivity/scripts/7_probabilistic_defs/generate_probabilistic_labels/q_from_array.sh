#!/bin/bash

selectionfile=/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/fmri_cohort/relmatch_parlabels3/namelist_full.txt

echo ${SGE_TASK_ID}
echo '*'

SUB=$(sed -n "$SGE_TASK_ID"p $selectionfile | awk '{print $1}')

echo $SUB
echo $RUNSTR

cd /home/weiner/shakki/scripts-github/analysis/misc/generate_probabilistic_labels

python generate_probabilistic_labels.py $SUB nat
python remove_label_overlap_mpm.py $SUB prob

