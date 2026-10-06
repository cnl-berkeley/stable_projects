### Assumes prediction_subs have been projected onto fsaverage and stored in
### $PREDICTION_DIR/projected_labels/$LABEL

SUBJECT_LIST=/home/weiner/HCP/subject_lists/HCP_processed_subs_all.txt

SUBJECTS=$(cat $SUBJECT_LIST)

SUBJECTS_DIR=/home/weiner/HCP/subjects/

PROJECT_DIR=/home/weiner/HCP/projects/ifrms_HCP/

PROJECT_ID=all_72_HCP

PREDICTION_SUBS=/home/weiner/HCP/subject_lists/HCP_processed_subs_all.txt

PREDICTION_DIR=/home/weiner/HCP/projects/ifrms_HCP/

LABELS='2'


## Create probability map directory if it doesn't exist already
PROB_MAP_DIR=${PROJECT_DIR}prob_maps
if [[ ! -e ${PROB_MAP_DIR} ]]; then
    mkdir ${PROB_MAP_DIR}
elif [[ ! -d ${PROB_MAP_DIR} ]]; then
    echo "$PROB_MAP_DIR already exists but is not a directory" 1>&2
fi


## Activate virtual env with relevant imports
echo ${SUB_ID}
cd /home/weiner/bparker/ifrms_HCP;
source activate /home/weiner/bparker/anaconda3/envs/python3_analysis_plotting;
unset PYTHONPATH;

for SUB_ID in $SUBJECTS; do

    
## Create subject prob_maps directory
src_label_dir=${PROJECT_DIR}prob_maps/${SUB_ID};
if [[ ! -e ${src_label_dir} ]]; then
       mkdir $src_label_dir; fi
    
## Python script creates MPM maps and MPM binary maps
# saved as ${PROB_MAP_DIR}/${SUB_ID}/${PROJECT_ID}_PROB_{${LABEL}.label

python /home/weiner/bparker/ifrms_HCP/MPM_labels.py ${SUB_ID} ${SUBJECTS_DIR} ${PROJECT_DIR} ${SUBJECT_LIST} ${PROJECT_ID} ${PREDICTION_DIR} ${PREDICTION_SUBS};

# register leave-one-out prob maps to left out subject 
SRC_SUB='fsaverage'   
TRG_SUBS=${SUB_ID}  



for TRG_SUB in ${TRG_SUBS};
do
  # make target directories 
  trg_label_dir=${SUBJECTS_DIR}${TRG_SUB}/label/fsaverage_MPM;
    
  if [[ ! -e ${trg_label_dir} ]]; then
      mkdir ${trg_label_dir};
  fi
  
 
  # BP: Commented out to remove merge labels; only interested in individual labels

  #mri_mergelabels -i ${src_label_dir}/lh.PROB_MPM_pmfs_p.label -i ${src_label_dir}/lh.PROB_MPM_pmfs_i.label -i ${src_label_dir}/lh.PROB_MPM_pmfs_a.label -o ${src_label_dir}/lh.PROB_MPM_mfs.label
  #mri_mergelabels -i ${src_label_dir}/rh.PROB_MPM_pmfs_p.label -i ${src_label_dir}/rh.PROB_MPM_pmfs_i.label -i ${src_label_dir}/rh.PROB_MPM_pmfs_a.label -o ${src_label_dir}/rh.PROB_MPM_mfs.label
  #for label in $labels

  
  for LABEL in $LABELS; do
    for filename in "${src_label_dir}/lh.${PROJECT_ID}_PROB_MPM_${LABEL}.label";
    do
     echo $filename
      
      trg_label_file=$(basename $filename)
      
      mri_label2label --srcsubject ${SRC_SUB} --srclabel ${filename} --trgsubject ${TRG_SUB} \
      --trglabel ${trg_label_dir}/${trg_label_file} --regmethod surface --hemi lh;
      
    done

    echo "lh MPMs transformed to ${TRG_sub} cortical surface";

    for filename in "${src_label_dir}/rh.${PROJECT_ID}_PROB_MPM_${LABEL}.label";
    do
      echo $filename
      
      trg_label_file=$(basename $filename)
      
      mri_label2label --srcsubject ${SRC_SUB} --srclabel ${filename} --trgsubject ${TRG_SUB} \
      --trglabel ${trg_label_dir}/${trg_label_file} --regmethod surface --hemi rh;
      
    done

    echo "rh MPMs transformed to ${TRG_sub} cortical surface";
    echo $trg_label_file
    
    
done
done

END_TIME=$(date);
echo "run completed successfully at $END_TIME"
done



