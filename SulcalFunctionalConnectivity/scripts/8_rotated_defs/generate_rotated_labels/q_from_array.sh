#!/bin/bash

selectionfile=/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/fmri_cohort/relmatch_parlabels3/namelist_full.txt

echo ${SGE_TASK_ID}
echo '*'

SUB=$(sed -n "$SGE_TASK_ID"p $selectionfile | awk '{print $1}')

echo $SUB
echo $RUNSTR

cd /home/weiner/shakki/scripts-new/null_models

python combine_labels_annot.py $SUB nat
python do_spins.py $SUB

