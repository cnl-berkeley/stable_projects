import os,sys
import os,sys
import nibabel as nib
from numpy import *
import gdist
#import surfdist as sd
import numpy as np


# NOTE: had trouble with numba format and jit in surfdist, so some functions are copied over with slight modifications below


def translate_src(src, cortex):
    """
    Convert source nodes to new surface (without medial wall).
    """
    src_new = np.array(np.where(np.in1d(cortex, src))[0], dtype=np.int32)

    return src_new


def triangles_keep_cortex(triangles, cortex):
    """
    Remove triangles with nodes not contained in the cortex label array
    """

    # for or each face/triangle keep only those that only contain nodes within the list of cortex nodes
    input_shape = triangles.shape
    triangle_is_in_cortex = np.all(np.reshape(np.in1d(triangles.ravel(), cortex), input_shape), axis=1)

    cortex_triangles_old = np.array(triangles[triangle_is_in_cortex], dtype=np.int32)

    # reassign node index before outputting triangles
    new_index = np.digitize(cortex_triangles_old.ravel(), cortex, right=True)
    cortex_triangles = np.array(np.arange(len(cortex))[new_index].reshape(cortex_triangles_old.shape), dtype=np.int32)

    return cortex_triangles


def surf_keep_cortex(surf, cortex):
    # split surface into vertices and triangles
    vertices, triangles = surf

    # keep only the vertices within the cortex label
    cortex_vertices = array(vertices[cortex], dtype=float64)

    # keep only the triangles within the cortex label
    cortex_triangles = triangles_keep_cortex(triangles, cortex)

    return cortex_vertices, cortex_triangles


def dist_calc_matrix(surf, cortex, label_inds_all):
    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)
    
    n_labels = len(labels)
    dist_mat = zeros((n_labels,n_labels))
    for r1 in arange(n_labels):
        for r2 in arange(n_labels):
            if any(label_inds_all[r1]) and any(label_inds_all[r2]):
                val2 = gdist.compute_gdist(cortex_vertices, cortex_triangles,
                                            source_indices = array(label_inds_all[r1]),
                                            target_indices = array(label_inds_all[r2]))
            dist_mat[r1,r2] = amin(val2)

    dist_mat = (dist_mat + dist_mat.T)/2

    return dist_mat


def getLabelIndices(sub,hemi,labels,cortex):
    label_inds_all = []

    n_labels = len(labels)
    for l in arange(n_labels):
        if type(labels[l]) is list: # pick the first label in list that exits
            label_found = False
            for lab in labels[l]:
                labelfile = '%s/sub-%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,lab)
                if os.path.exists(labelfile) and not label_found:
                    labelfile_use = labelfile
                    label_found = True
                else:
                    labelfile = '%s/sub-%s/sulcal_prob_mpm/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,lab)
                    if os.path.exists(labelfile) and not label_found:
                        labelfile_use = labelfile
                        label_found = True

        else: # look for specific label
            label_found = False
            labelfile = '%s/sub-%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,labels[l])
            if os.path.exists(labelfile) and not label_found:
                labelfile_use = labelfile
                label_found = True
            else:
                labelfile = '%s/sub-%s/sulcal_prob_mpm/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,labels[l])
                if os.path.exists(labelfile) and not label_found:
                    labelfile_use = labelfile
                    label_found = True

            labelfile_use = '%s/sub-%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,labels[l])

        label_inds = nib.freesurfer.io.read_label(labelfile_use, read_scalars=False)
        label_inds_t = translate_src(label_inds,cortex) # exclude medial wall
        label_inds_all.append(label_inds_t)
    
    return label_inds_all


def getDistMatrix(fmriprep_derivatives_dir,annot_file,sub,hemi,savedir):
    highres_surface = '%s/sub-%s/ses-%s/anat/sub-%s_ses-%s_hemi-%s_midthickness.surf.gii'%(fmriprep_derivatives_dir,sub,sub[-1],sub,sub[-1],hemi[0].upper())
    giidata = nib.load(highres_surface)
    giidata2 = squeeze(asarray([x.data for x in giidata.darrays])) 
    surf = (giidata2[0],giidata2[1])  

    cort_file = '%s/sub-%s/label/%s.cortex.label'%(os.environ['SUBJECTS_DIR'],sub,hemi)
    cortex = sort(nib.freesurfer.read_label(cort_file))

    label_inds_all = getLabelIndices(sub,hemi,labels,cortex)

    dist_matrix = dist_calc_matrix(surf,cortex,label_inds_all)
    print('dist_matrix',dist_matrix)

    savetxt('%s/adj-labels-%s.txt'%(savedir,hemi),dist_matrix)


if __name__ == '__main__':

    sub = sys.argv[1]
    hemi = sys.argv[2]
    prob_tag = sys.argv[3]

    if prob_tag != 'prob':
        outdir='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/labelnw'%sub
        mpar_labels = []
        lpfc_labels = ['ifs',['painfs_any','painfs_combined'],'pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
        lpar_labels = [['slos1','slocs-v','SLOS'],'sB','pips','mTOS',['iTOS','ITOS','lTOS'],'IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
        labels = lpfc_labels + lpar_labels + mpar_labels
    else:
        outdir='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/labelnw_prob_mpm'%sub
        lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
        lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
        mpar_labels = []
        labels = lpfc_labels + lpar_labels + mpar_labels
        for l in np.arange(len(labels)):
              labels[l] = labels[l] + '_prob33'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    os.environ['SUBJECTS_DIR'] = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7'
    fmriprep_derivatives_dir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/bids/derivatives_v7'

    getDistMatrix(fmriprep_derivatives_dir,labels,sub,hemi,outdir)
