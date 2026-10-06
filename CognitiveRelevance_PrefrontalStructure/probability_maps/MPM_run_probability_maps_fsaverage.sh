
# Path to subjects list
 SUBJECT_LIST=/home/weiner/HCP/subject_lists/HCP_processed_subs_all.txt

# Subjects
 SUBJECTS=$(cat $SUBJECT_LIST)

# Subjects directory
 SUBJECTS_DIR=/home/weiner/HCP/subjects

# Projection directory
 PROJECT_DIR=/home/weiner/HCP/subjects/fsaverage/painfs

# Name of project
 PROJECT_ID=PIMFS

# Participants
 PREDICTION_SUBS=/home/weiner/HCP/subject_lists/HCP_processed_subs_all.txt

# Fsaverage
 PREDICTION_DIR=/home/weiner/HCP/subjects/fsaverage/painfs

# Sulcal labels
 LABELS[1]='pimfs_d'
 LABELS[2]='pimfs_v'




## Create probability map directory if it doesn't exist already
#PROB_MAP_DIR=${PROJECT_DIR}/prob_maps
#if [[ ! -e ${PROB_MAP_DIR} ]]; then
#mkdir ${PROB_MAP_DIR}
#elif [[ ! -d ${PROB_MAP_DIR} ]]; then
    #echo "$PROB_MAP_DIR already exists but is not a directory" 1>&2
#fi

## Activate virtual env with relevant imports
#echo ${SUB_ID}
#cd /home/weiner/bparker/ifrms_HCP;
#source activate /home/weiner/bparker/anaconda3/envs/python3_analysis_plotting;
#unset PYTHONPATH;

SUB_ID=fsaverage
echo $SUB_ID

    
## Create subject prob_maps directory
#src_label_dir=${PROJECT_DIR}/prob_maps/${SUB_ID};
#if [[ ! -e ${src_label_dir} ]]; then
       #mkdir $src_label_dir; fi
    
## Python script creates MPM maps and MPM binary maps
## saved as ${PROB_MAP_DIR}/${SUB_ID}/${PROJECT_ID}_PROB_{${LABEL}.label

## NOTE: IN MPM_label.py UPDATE ALL LABEL NAMES


#for SUB_ID in $SUBJECTS; do
#SUB_ID='fsaverage'
    
## Create subject prob_maps directory
#src_label_dir=${PROJECT_DIR}/prob_maps/${SUB_ID};
#if [[ ! -e ${src_label_dir} ]]; then
      # mkdir $src_label_dir; fi
    
## Python script creates MPM maps and MPM binary maps
## saved as ${PROB_MAP_DIR}/${SUB_ID}/${PROJECT_ID}_PROB_{${LABEL}.label

## NOTE: IN MPM_label.py UPDATE ALL LABEL NAMES

python /home/weiner/HCP/projects/HCP_reasoning/pimfs/mpm/MPM_labels_fsaverage.py ${SUB_ID} ${SUBJECTS_DIR} ${PROJECT_DIR} ${SUBJECT_LIST} ${PROJECT_ID} ${PREDICTION_DIR} ${PREDICTION_SUBS};


END_TIME=$(date);
echo "run completed successfully at $END_TIME"
#done
