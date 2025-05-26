#!/bin/bash

selectionfile=/home/weiner/shakki/scripts-new/analysis/group-analysis/ml-prob/sklearn_mpm/run_analyses.txt

# Resolve number of inputs, according to the lines in namelist.txt file
no_folders=(`wc -l $selectionfile`)
no_folder=${no_folders[0]}
echo $no_folders " analyses listed"

qsub -pe threaded 3 -binding linear:3 -N classification -t 1-${no_folders} -j y ./q_from_array.sh

