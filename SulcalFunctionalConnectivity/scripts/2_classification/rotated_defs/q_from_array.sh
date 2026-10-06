#!/bin/bash

echo ${SGE_TASK_ID}

cd /home/weiner/shakki/scripts-github/analysis/classification/chance_level_analysis

python do-svm-pair-scale-log5-perm.py $SGE_TASK_ID

