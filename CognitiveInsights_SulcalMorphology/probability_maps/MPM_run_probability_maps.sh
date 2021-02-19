## This script creates maximum probability maps for a given subject based labels by the same name from the subjects in PREDICTION_SUBS ###
## NOTE: if PREDICTION_SUBS and SUBJECTS are the same, the maps will leave out one subject and predict based on all others within SUBJECTS ###

## Assumes all individual subject labels have been projected to fsaverage ##




SUBJECT_LIST=###SUBJECT LIST .TXT FOR THE TARGET SUBJECTS OF PREDICTIONS

SUBJECTS=$(cat $SUBJECT_LIST)

SUBJECTS_DIR=### FREESURFER SUBJECT DIRECTORY FOR SUBJECTS IN SUBJECT_LIST

PROJECT_DIR=### DIRECTORY WITH FSAVERAGE PROJECTED LABELS FOR TARGET SUBJECTS

PROJECT_ID=### PROJECT ID, WILL BE ADDED TO OUTPUT FILE

PREDICTION_SUBS=### SUBJECT LIST FOR SOURCE OF PREDICTIONS (New subject list if predicting on second dataset. Otherwise same as SUBJECTS_LIST)

PREDICTION_DIR=### DIRECTORY WITH FSAVERAGE PROJECTED LABELS FOR PREDICTION SUBJECTS (New directory if predicting on second dataset. Otherwise same as PROJECT_DIR)

LABELS=### LIST OF SULCAL LABEL NAMES


## Create probability map directory if it doesn't exist already
PROB_MAP_DIR=${PROJECT_DIR}prob_maps
if [[ ! -e ${PROB_MAP_DIR} ]]; then
    mkdir ${PROB_MAP_DIR}
elif [[ ! -d ${PROB_MAP_DIR} ]]; then
    echo "$PROB_MAP_DIR already exists but is not a directory" 1>&2
fi


for SUB_ID in $SUBJECTS; do

    
## Create subject prob_maps directory
src_label_dir=${PROJECT_DIR}prob_maps/${SUB_ID};
if [[ ! -e ${src_label_dir} ]]; then
       mkdir $src_label_dir; fi
    
## Python script creates MPM maps and MPM binary maps
## saved as ${PROB_MAP_DIR}/${SUB_ID}/${PROJECT_ID}_PROB_{${LABEL}.label

## NOTE: IN MPM_label.py UPDATE ALL LABEL NAMES

python /home/weiner/bparker/ifrms_HCP/MPM_labels.py ${SUB_ID} ${SUBJECTS_DIR} ${PROJECT_DIR} ${SUBJECT_LIST} ${PROJECT_ID} ${PREDICTION_DIR} ${PREDICTION_SUBS} ${LABELS};

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



