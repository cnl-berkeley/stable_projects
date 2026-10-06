#!/bin/bash

DERIVATIVES_DIR=/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/bids/derivatives_v7
BIDS_DIR=/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/bids
FS_SUBJECTS_DIR=/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7

OUTDIR=/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/postproc_fmri_smooth0mm_v7-bp/relmatch

selectionfile=/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/postproc_fmri/relmatch-bp/namelist.txt

#echo $selectionfile
echo ${SGE_TASK_ID}
echo '*'

SUB=$(sed -n "$SGE_TASK_ID"p $selectionfile | awk '{print $1}')
RUNSTR="0,1,2"

echo $SUB
echo $RUNSTR
WORK_DIR=/home/weiner/shakki/NotBackedUp/postproc_fmri/relmatch/tmp-relmatch-${SGE_TASK_ID}-${SUB}

cd /home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/postproc_fmri/relmatch-bp
python run-denoise.py $DERIVATIVES_DIR $OUTDIR --participant-label $SUB --runstr $RUNSTR --task relmatch --bids-dir $BIDS_DIR --fs_subjects_dir $FS_SUBJECTS_DIR --work-dir $WORK_DIR --smooth_surf_mm 0 --acompcor --friston24

rm -rf $WORK_DIR
