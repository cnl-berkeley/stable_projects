import os, sys

import numpy as np

from netneurotools import freesurfer as nnsurf

import nibabel as nib 


def get_coords_native(sub):
    lhsphere = '%s/%s/surf/lh.sphere.reg'%(os.environ['SUBJECTS_DIR'],sub)
    rhsphere = '%s/%s/surf/rh.sphere.reg'%(os.environ['SUBJECTS_DIR'],sub)

    coords, hemi = [], []
    for n, sphere in enumerate([lhsphere, rhsphere]):
        coords.append(nib.freesurfer.read_geometry(sphere)[0])
        hemi.append(np.ones(len(coords[-1])) * n)

    return np.vstack(coords), np.hstack(hemi)


def projectProbMask(prob_label,sub,hemi):
    from nipype.interfaces.freesurfer import Label2Label
    l2l = Label2Label()
    l2l.inputs.hemisphere = hemi
    l2l.inputs.subject_id = sub
    l2l.inputs.sphere_reg = '%s/%s/surf/%s.pial'%(os.environ['SUBJECTS_DIR'],sub,hemi)
    l2l.inputs.white = '%s/%s/surf/%s.pial'%(os.environ['SUBJECTS_DIR'],sub,hemi)
    l2l.inputs.source_subject = 'fsaverage_prob'
    l2l.inputs.source_label = prob_label
    l2l.inputs.source_white = '%s/%s/surf/%s.pial'%(os.environ['SUBJECTS_DIR'],'fsaverage_prob',hemi)
    l2l.inputs.source_sphere_reg = '%s/%s/surf/%s.pial'%(os.environ['SUBJECTS_DIR'],'fsaverage_prob',hemi)
    stem = prob_label.split('/')[-1].split('.')[1]
    l2l.inputs.out_file = '%s/%s/label/%s.%s_converted.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,stem)
    res = l2l.run()
    print(res.outputs)

    return res.outputs.out_file


def _decode_list(vals):
    """List decoder."""
    return [val.decode() if hasattr(val, 'decode') else val for val in vals]


try:  # scipy >= 1.8.0
    from scipy.ndimage._measurements import _stats, labeled_comprehension
except ImportError:  # scipy < 1.8.0
    from scipy.ndimage.measurements import _stats, labeled_comprehension

def spin_parcels_native(spins, lhannot, rhannot, lhmask, rhmask, n_rotate=2, drop=None):
    def overlap(vals):
        """Return most common non-negative value in `vals`; -1 if all neg."""
        vals = np.asarray(vals)
        vals, counts = np.unique(vals[vals > 0], return_counts=True)
        try:
            return vals[counts.argmax()]
        except ValueError:
            return -1

    regional_inds_lh = nib.freesurfer.io.read_label(lhmask)
    regional_inds_rh = nib.freesurfer.io.read_label(rhmask)
    regional_inds = [regional_inds_lh,regional_inds_rh]
    drop = ['lh.not_in_region','rh.not_in_region']

    # get vertex-level labels (set drop labels to - values)
    vertices, end = [], 0
    for n, annot in enumerate([lhannot, rhannot]):
        # modify annotation
        if n == 0:
            hemi = 'lh'
        else:
            hemi = 'rh'
        labels, ctab, names = nib.freesurfer.read_annot(annot)
        vertices_hemi = np.arange(len(labels))
        names2 = []
        for i in np.arange(len(names)):
            names2.append('%s.%s'%(hemi,names[i].decode('UTF-8')))
        names = list(names2)
        labels[~np.isin(vertices_hemi,regional_inds[n])] = np.amax(labels)+1 # to be omitted
        labels = labels + len(np.unique(vertices))
        names.append('%s.not_in_region'%hemi)

        # combine info on hemis and regions to be omitted
        names = _decode_list(names)
        todrop = set(names) & set(drop)
        inds = [names.index(f) - n for n, f in enumerate(todrop)]
        labs = np.arange(len(names) - len(inds)) + (end - (len(inds) * n))
        insert = np.arange(-1, -(len(inds) + 1), -1)
        vertices.append(np.insert(labs, inds, insert)[labels])
        end += len(names)
    vertices = np.hstack(vertices)
    labels = np.unique(vertices)
    mask = labels > -1

    # spin and assign regions based on max overlap
    regions = np.zeros((len(labels[mask]), n_rotate), dtype='int32')
    for n in range(n_rotate):
        msg = f'Calculating parcel overlap: {n:>5}/{n_rotate}'
        print(msg, end='\b' * len(msg), flush=True)
        regions[:, n] = labeled_comprehension(vertices[spins[:, n]], vertices,
                                              labels, overlap, int, -1)[mask]

    print(regions.shape)

    return regions


def saveLabel(verts_sel,coords,savepath):
    of = open(savepath,'w')
    n_label_verts = len(verts_sel)
    of.write('%d\n'%n_label_verts)
    for v in np.arange(n_label_verts):
        vert = verts_sel[v]
        of.write('%d  %.3f  %.3f  %.3f  %.6f\n'%(vert,coords[vert,0],coords[vert,1],coords[vert,2],0.0))
    of.write('\n')
    of.close()


if __name__ == '__main__':
    # modified to native surfaces from
    # https://markello-spatialnulls.netlify.app/spatial_nulls.html

    os.environ['SUBJECTS_DIR'] = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7'
    sub = 'sub-%s'%sys.argv[1]

    outdir = '%s/%s/rotated_labels'%(os.getcwd(),sub)
    if not os.path.exists(outdir):
         os.makedirs(outdir)

    # load surfaces
    coords, hemi = get_coords_native(sub)

    # load sulcal definitions
    lhannot = '%s/%s/label/lh.LPFC+LPC.annot'%(os.environ['SUBJECTS_DIR'],sub)
    rhannot = '%s/%s/label/rh.LPFC+LPC.annot'%(os.environ['SUBJECTS_DIR'],sub)
    lhlabels, ctab, lhnames = nib.freesurfer.read_annot(lhannot)
    rhlabels, ctab, rhnames = nib.freesurfer.read_annot(rhannot)
    n_lh_labels = np.max(lhlabels)
    n_lh_vert = len(lhlabels)
    rhlabels[rhlabels > 0] += n_lh_labels
    labels = np.hstack([lhlabels,rhlabels])

    uni_labels = np.unique(labels[labels > 0])
    n_labels = len(uni_labels)
    print('n_labels',n_labels)
    print('labels',labels)
    
    labelnames = []
    for i in np.arange(len(lhnames)):
        labelnames.append('lh.%s'%lhnames[i].decode('UTF-8'))
    for i in np.arange(len(rhnames)):
        labelnames.append('rh.%s'%rhnames[i].decode('UTF-8'))
    labelnames = np.array(labelnames)
    labelnames = labelnames[labelnames != 'lh.Unknown']
    labelnames = labelnames[labelnames != 'rh.Unknown']
    labeltypes = ['LPFC']*10 + ['LPC']*11 + ['LPFC']*10 + ['LPC']*11
    labelhemi = ['lh']*21 + ['rh']*21
    print(labelnames)

    # per label, filter to rotations within LPFC
    prob_LPFC = ['/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7/fsaverage_prob/label/lh.LPFC.label',\
    '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7/fsaverage_prob/label/rh.LPFC.label']
    lh_LPFC_mask = projectProbMask(prob_LPFC[0],sub,'lh')
    rh_LPFC_mask = projectProbMask(prob_LPFC[1],sub,'rh')
#    LPFC_checks = spin_parcels_native(spins, lhannot, rhannot, lh_LPFC_mask, rh_LPFC_mask)
#    print(LPFC_checks.shape)
 
    prob_LPC = ['/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7/fsaverage_prob/label/lh.LPC.label',\
    '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7/fsaverage_prob/label/rh.LPC.label']
    lh_LPC_mask = projectProbMask(prob_LPC[0],sub,'lh')
    rh_LPC_mask = projectProbMask(prob_LPC[1],sub,'rh')

    n_vertices = coords.shape[0] #np.arange(coords.shape[0])

    lpfc_inds_lh = nib.freesurfer.io.read_label(lh_LPFC_mask)
    lpfc_inds_rh = nib.freesurfer.io.read_label(rh_LPFC_mask)
    lpfc_inds_rh += n_lh_vert
    lpfc_inds = np.hstack([lpfc_inds_lh,lpfc_inds_rh])
    lpfc_test_data = np.ones_like(labels)*np.nan
    lpfc_test_data[lpfc_inds] = 1

    lpc_inds_lh = nib.freesurfer.io.read_label(lh_LPC_mask)
    lpc_inds_rh = nib.freesurfer.io.read_label(rh_LPC_mask)
    lpc_inds_rh += n_lh_vert
    lpc_inds = np.hstack([lpc_inds_lh,lpc_inds_rh])
    lpc_test_data = np.ones_like(labels)*np.nan
    lpc_test_data[lpc_inds] = 1

    print(np.unique(lpfc_inds))
    print(np.unique(lpc_inds))

    ready = False
    n_wanted_samples = 1000
    c = np.zeros((n_labels,1))
    print(c)

    rot_dict = {}
    for i in np.arange(n_labels):
        rot_dict[i] = []
    counter = 0
    while not ready:
        # generate rotations
        n_rotate = 500
        spins = nnsurf.gen_spinsamples(coords, hemi, exact=False,
                                       n_rotate=n_rotate, verbose=True,
                                       seed=counter, check_duplicates=False)

        # apply rotations
        annot_rots = np.zeros((len(labels),n_rotate), dtype='int32')
        for n in np.arange(n_rotate):
            annot_rots[:,n] = labels[spins[:,n]]

        # check rotations
        for n in np.arange(n_rotate):
            for l in np.arange(n_labels):
                rot_label_inds = np.argwhere(annot_rots[:,n] == uni_labels[l])
                rot_label_inds = [j for i in rot_label_inds for j in i]
                if labeltypes[l] == 'LPFC':
                     test_data = lpfc_test_data
                else:
                     test_data = lpc_test_data

                if not np.isnan(np.mean(test_data[rot_label_inds])):
                     rot_dict[l].append(rot_label_inds)
                     c[l] += 1

        if np.min(c,0) >= n_wanted_samples:
            ready = True

        counter+= 1

        print(counter, c.T)

#    np.savez('rot_dict.npz',rot_dict=rot_dict)
#    npz = np.load('rot_dict.npz',allow_pickle=True)
#    rot_dict = npz['rot_dict']
    for l in np.arange(n_labels):
       of = open('%s/%s_rotations.csv'%(outdir,labelnames[l]),'w')
       for i in np.arange(n_wanted_samples):
           for j in rot_dict[l][i]:
               of.write('%s,'%j)
           of.write('\n')
       of.close()

    """
    for l in np.arange(n_labels):
       #arr = rot_dict[l]
       arr = np.loadtxt('%s/%s_rotations.csv'%(outdir,labelnames[l])
       # save test overlay with the number of hits per vertex
       vec = np.zero((n_vertices,1))
       vals_flat = [x for xs in arr for x in xs]
       values,counts = np.unique(rot_dict[l],return_counts=True)
       for v in np.arange(len(values)):
           vec[values[v]] = counts[v]
       tmp = vec.astype('>f4')

       if labelhemi[l] == 0:
           mgh_lh = nib.freesurfer.mghformat.MGHImage(vec[0:n_lh_vert,0].astype('>f4'),np.eye(4))
           mgh_lh.to_filename('%s/%srothits_lh.mgh'%(outdir,labelnames[l]))
       else:
           vec = vec[n_lh_vert:len(labels)]
           vec = vec - n_lh_vert
           mgh_rh = nib.freesurfer.mghformat.MGHImage(vec.astype('>f4'),np.eye(4))
           mgh_rh.to_filename('%s/%s_rothits_rh.mgh'%(outdir,labelnames[l]))
   """
