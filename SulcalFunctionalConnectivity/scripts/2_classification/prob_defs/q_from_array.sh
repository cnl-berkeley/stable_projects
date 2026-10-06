#!/bin/bash

selectionfile=/home/weiner/shakki/scripts-new/analysis/group-analysis/ml-prob/sklearn_mpm/run_analyses.txt

echo $selectionfile
echo ${SGE_TASK_ID}
echo '*'

NPZ=$(sed -n "$SGE_TASK_ID"p $selectionfile | awk '{print $1}')
SAVENAME=$(sed -n "$SGE_TASK_ID"p $selectionfile | awk '{print $2}')

echo $NPZ
echo $SAVENAME

cd /home/weiner/shakki/scripts-new/analysis/group-analysis/ml-prob/sklearn_mpm

python do-svm-pair-scale-log5-perm.py $NPZ $SAVENAME

