# Set subject directory (assumes fsaverage surface is located there)
SUBJECTS_DIR=

# Choose subject(s) to project label to
TRG_SUBS=''

# Set label directory to location of probability maps
PROB_MAP_DIR=

# Select which probability map you want to project to individual-level surfaces (e.g., thresholded: 'PCGS_Urgency_lh_PROB_MPM_binary_0.25_PCGS' or unthresholded: 'PCGS_Urgency_lh_PROB_PCGS')
PROB_MAP=''

# Run mri_label2label to project to new surface
for TRG_SUB in $TRG_SUBS;
do

## Left hemisphere
mri_label2label --srcsubject $SUBJECTS_DIR/fsaverage --srclabel $PROB_MAP_DIR/lh.$PROB_MAP=''.label --trgsubject $SUBJECTS_DIR/$TRG_SUB --trglabel $SUBJECTS_DIR/$TRG_SUB/label/lh.$PROB_MAP=''.label --regmethod surface

## Right hemisphere
mri_label2label --srcsubject $SUBJECTS_DIR/fsaverage --srclabel $PROB_MAP_DIR/rh.$PROB_MAP=''.label --trgsubject $SUBJECTS_DIR/$TRG_SUB --trglabel $SUBJECTS_DIR/$TRG_SUB/label/rh.$PROB_MAP=''.label --regmethod surface

done