#!/bin/bash

f=$(sed -n "$SGE_TASK_ID"p $RUNLIST | awk '{print $1}')

mkdir -p $LOGDIR

fmriprep-docker ${BIDS_DIR} ${BIDS_DIR}/derivatives participant --participant-label $f --fs-license-file ${FS_LIC_PATH}/license.txt  --output-spaces func anat fsnative MNI152NLin2009cAsym --dummy-scans 3 --skull-strip-t1w force --fs-subjects-dir $FS_SUBJECTS_DIR --notrack -vv --skip_bids_validation
