## Anchoring functional connectivity to individual sulcal morphology yields insights in a pediatric study of reasoning ##

Analysis code associated with Häkkinen et al. (2025) *Anchoring functional connectivity to individual sulcal morphology yields insights in a pediatric study of reasoning*. The published manuscript is available at The Journal of Neuroscience.

For additional questions or data requests please email Silvia Bunge (sbunge@berkeley.edu).

### Scripts ###

Code relating to different subtasks is organized to separate folders: connectivity estimation, classification, clustering, network correlations, graph metrics, depth correlations, and miscellaneous. Scripts are in Python and Matlab, and rely on various open source toolboxes (e.g. FreeSurfer, NiPype, scikit-learn); the necessary libaries are listed in the beginning of each script.

Input is expected to include fMRIprep-preprocessed functional data, cortical reconstructions by FreeSurfer, and sulcal definitions in FreeSurfer label format.

### Associated data ###

Probability maps: Spatial probability maps for each of 42 LPFC-LPC sulci can be found in the folder probability_maps. The maps are in standard FreeSurfer fsaverage space. Values are unthresholded, showing the ratio of subjects (N = 43) labeled for the sulcus per vertex.
