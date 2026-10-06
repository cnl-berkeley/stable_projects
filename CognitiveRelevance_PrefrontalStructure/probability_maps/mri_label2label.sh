# Set subject directory (assumes fsaverage surface is located there)
SUBJECTS_DIR=

# Choose subject(s) to project label to
TRG_SUBS=''

# Set label directory to location of pimfs probability maps
PROB_MAP_DIR=

# Rename if you want to project a different probability map from folder (e.g., threshold = 0.2 or unthresholded)
PROB_MAP_PIMFS_D='PIMFS_PROB_MPM_0.1_pimfs_d'
PROB_MAP_PIMFS_V='PIMFS_PROB_MPM_0.1_pimfs_v'

# Run mri_label2label to project to new surface
for TRG_SUB in $TRG_SUBS;
do

## Pimfs-d left hemisphere
mri_label2label --srcsubject $SUBJECTS_DIR/fsaverage --srclabel $PROB_MAP_DIR/lh.$PROB_MAP_PIMFS_D.label --trgsubject $SUBJECTS_DIR/$TRG_SUB --trglabel $SUBJECTS_DIR/$TRG_SUB/label/lh.$PROB_MAP_PIMFS_D.label --regmethod surface

## Pimfs-v left hemisphere
mri_label2label --srcsubject $SUBJECTS_DIR/fsaverage --srclabel $PROB_MAP_DIR/lh.$PROB_MAP_PIMFS_V.label --trgsubject $SUBJECTS_DIR/$TRG_SUB --trglabel $SUBJECTS_DIR/$TRG_SUB/label/lh.$PROB_MAP_PIMFS_V.label --regmethod surface

## Pimfs-d right hemisphere
mri_label2label --srcsubject $SUBJECTS_DIR/fsaverage --srclabel $PROB_MAP_DIR/rh.$PROB_MAP_PIMFS_D.label --trgsubject $SUBJECTS_DIR/$TRG_SUB --trglabel $SUBJECTS_DIR/$TRG_SUB/label/rh.$PROB_MAP_PIMFS_D.label --regmethod surface

## Pimfs-v right hemisphere
mri_label2label --srcsubject $SUBJECTS_DIR/fsaverage --srclabel $PROB_MAP_DIR/rh.$PROB_MAP_PIMFS_V.label --trgsubject $SUBJECTS_DIR/$TRG_SUB --trglabel $SUBJECTS_DIR/$TRG_SUB/label/rh.$PROB_MAP_PIMFS_V.label --regmethod surface

done