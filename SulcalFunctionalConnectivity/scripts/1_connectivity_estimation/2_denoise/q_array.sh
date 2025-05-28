#!/bin/bash

selectionfile=/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/postproc_fmri/relmatch-bp/namelist.txt

# Resolve number of inputs, according to the lines in namelist.txt file
no_folders=(`wc -l $selectionfile`)
no_folder=${no_folders[0]}
echo $no_folders " subjects to be processed"

qsub -pe threaded 3 -binding linear:3 -N postproc-rsfmri -t 1-${no_folders} -j y ./q_from_array.sh

