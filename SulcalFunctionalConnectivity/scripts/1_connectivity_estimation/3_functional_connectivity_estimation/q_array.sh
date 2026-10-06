#!/bin/bash

selectionfile=/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/fmri_cohort/relmatch_parlabels3/namelist_full.txt

# Resolve number of inputs, according to the lines in namelist.txt file
no_folders=(`wc -l $selectionfile`)
no_folder=${no_folders[0]}
echo $no_folders " subjects to be processed"

MEM=2G
qsub -pe threaded 3 -l h_vmem=$MEM -l mem_free=$MEM -binding linear:3 -N fc-est -t 1-${no_folders} -j y ./q_from_array.sh

