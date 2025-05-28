import os,sys
import nibabel as nib
import numpy as np
from copy import deepcopy


def registerLabel(label_file,labelname,hemi,tsub,ssub,outdir):
    from nipype.interfaces.freesurfer import Label2Label
    l2l = Label2Label()
    l2l.inputs.hemisphere = hemi
    l2l.inputs.subject_id = tsub #'sub-%s'%tsub
    l2l.inputs.sphere_reg = '%s/%s/surf/%s.sphere.reg'%(os.environ['SUBJECTS_DIR'],tsub,hemi)
    l2l.inputs.white = '%s/%s/surf/%s.white'%(os.environ['SUBJECTS_DIR'],tsub,hemi)
    l2l.inputs.source_sphere_reg = '%s/%s/surf/%s.sphere.reg'%(os.environ['SUBJECTS_DIR'],ssub,hemi)
    l2l.inputs.source_white = '%s/%s/surf/%s.white'%(os.environ['SUBJECTS_DIR'],ssub,hemi)
    l2l.inputs.source_subject = ssub #'sub-%s'%ssub
    l2l.inputs.source_label = label_file
    l2l.inputs.out_file = '%s/%s.%s_%s.label'%(outdir,hemi,labelname,ssub)
    res = l2l.run()

    return res.outputs.out_file


def getLabelIndsNat(labels,labelname,hemi,subs,sub,outdir):
    if sub != 'fsaverage':
        sub = 'sub-%s'%sub
    label_inds_nat = []
    for s in np.arange(len(subs)):
        reg_file = registerLabel(labels[s],labelname,hemi,sub,'sub-%s'%subs[s],outdir)
        label_inds_nat.append(list(nib.freesurfer.read_label(reg_file)))

    return label_inds_nat


def getLabelFile(sub,hemi,label):
    labelfile = ''
    file_found = ''
    if type(label) is list: # pick the first label in list that exits
        label_found = False
        for lab in label:
            labelfile = '%s/sub-%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,lab)
            if os.path.exists(labelfile) and not label_found:
                file_found = labelfile
                label_found = True
    else: # look for specific label
        labelfile = '%s/sub-%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,label)
        file_found = labelfile

    return file_found


def saveProbabilisticOverlay(label_inds_sel,coords,savepath):
    # https://surfer.nmr.mgh.harvard.edu/fswiki/LabelsClutsAnnotationFiles
    n_verts = len(coords)

    n_subs = len(label_inds_sel)
    
    arr = np.zeros((n_subs,n_verts))
    for i in np.arange(n_subs):
        arr[i,label_inds_sel[i]] = 1
    sumvec = np.sum(arr,axis=0)
    print(sumvec)
    prob_vec = sumvec/n_subs
    prob_vec = prob_vec.astype('>f4')

    mgh_img = nib.freesurfer.mghformat.MGHImage(prob_vec,np.eye(4))
    mgh_img.to_filename(savepath)

    return prob_vec


def saveProbLabelForAll(subs,hemi,label,outdir):
    n_subs = len(subs)
    
    labelname = ''
    if type(label) is list: # pick the first label in list that exits
        labelname = label[0]
    else:
        labelname = label

    # collect label files for all subjects
    label_files_all = []
    for s in np.arange(n_subs):
        label_files_all.append(getLabelFile(subs[s],hemi,label))
    print(label_files_all)

    # for each, compute probabilistic label based on all other subjects
    for s in np.arange(n_subs):
        print(subs[s])

        savedir = '%s/new/sub-%s'%(outdir,subs[s])
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        savedir_tmp = '/home/weiner/shakki/NotBackedUp/prob_label_tmp/tmp'
        if not os.path.exists(savedir_tmp):
            os.makedirs(savedir_tmp)

        label_files_sel = deepcopy(label_files_all)
        del label_files_sel[s] # do not include the subject

        subs_sel = deepcopy(subs)
        del subs_sel[s]

        # collect label indices in the native spece of this subject
        label_inds_sel = getLabelIndsNat(label_files_sel,labelname,hemi,subs_sel,subs[s],savedir_tmp)

        # coordinates to save the label file
        surf_file = '%s/sub-%s/surf/%s.midthickness'%(os.environ['SUBJECTS_DIR'],subs[s],hemi)
        coords,faces = nib.freesurfer.io.read_geometry(surf_file)

        # save probabilistic overlay
        prob_vec = saveProbabilisticOverlay(label_inds_sel,coords,'%s/%s_%s_prob.mgh'%(savedir,hemi,labelname))

        # save label at 33% probability
        prob_vec[prob_vec <= 0.33] = 0
        verts_sel = np.argwhere(prob_vec > 0) # vertex indices in label (0-based!)
        of = open('%s/%s.%s_prob33_mpm.label'%(savedir,hemi,labelname),'w') # to SUBJECTS_DIR?
        n_label_verts = len(verts_sel)
        of.write('%d\n'%n_label_verts)
        for v in np.arange(n_label_verts):
            vert = verts_sel[v]
            of.write('%d  %.3f  %.3f  %.3f  %.6f\n'%(vert,coords[vert,0],coords[vert,1],coords[vert,2],prob_vec[vert]))
        of.write('\n')
        of.close()


def saveProbLabelAverage(subs,hemi,label,outdir):
    n_subs = len(subs)
    
    labelname = ''
    if type(label) is list: # pick the first label in list that exits
        labelname = label[0]
    else:
        labelname = label

    # collect label files for all subjects
    label_files_all = []
    for s in np.arange(n_subs):
        label_files_all.append(getLabelFile(subs[s],hemi,label))
    print(label_files_all)

    # fsaverage - use all subjects
    sub = 'fsaverage'

    savedir = '%s/new/%s'%(outdir,sub)
    if not os.path.exists(savedir):
         os.makedirs(savedir)
    savedir_tmp = '%s/tmp'%savedir
    if not os.path.exists(savedir_tmp):
        os.makedirs(savedir_tmp)

    # collect label indices in the native spece of this subject
    label_inds_all = getLabelIndsNat(label_files_all,labelname,hemi,subs,sub,savedir_tmp)

    # coordinates to save the label file
    surf_file = '%s/%s/surf/%s.white'%(os.environ['SUBJECTS_DIR'],sub,hemi)
    coords,faces = nib.freesurfer.io.read_geometry(surf_file)

    # save probabilistic overlay
    prob_vec = saveProbabilisticOverlay(label_inds_all,coords,'%s/%s_%s_prob.mgh'%(savedir,hemi,labelname))

    # save label at 33% probability
    prob_vec[prob_vec <= 1.0/100] = 0
    verts_sel = np.argwhere(prob_vec > 0) # vertex indices in label (0-based!)
    of = open('%s/%s.%s_full_mpm.label'%(savedir,hemi,labelname),'w')
    of.write('#!ascii label, from subject fsaverage, vox2ras=TkReg\n')
    n_label_verts = len(verts_sel)
    of.write('%d\n'%n_label_verts)
    for v in np.arange(n_label_verts):
        vert = verts_sel[v]
        of.write('%d  %.3f  %.3f  %.3f %.6f\n'%(vert,coords[vert,0],coords[vert,1],coords[vert,2],prob_vec[vert]))
    of.write('\n')
    of.close()


if __name__ == '__main__':
    outdir = os.getcwd()

    lpfc_labels = ['ifs','painfs_a','painfs_p',['painfs_any','painfs_combined'],'pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
    lpar_labels = [['STS','sts'],['SPS','sps'],['SmgS','smgS','Smgs'],['slos1','slocs-v','SLOS'],'sB','pips','mTOS',['iTOS','ITOS','lTOS'],'IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    mpar_labels = []
    labels = lpfc_labels + lpar_labels + mpar_labels

    os.environ['SUBJECTS_DIR'] = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7'

    pipeline = 'orig'
    bin_tag = 'bin'
    hemis = 'both'
    pipe = '-bp'
    norm = 'nonorm'
    tag = '-sareg-fdreg'
    subs = np.genfromtxt('/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-%s/regs/namelist.csv'%(pipe,pipeline,tag,norm,hemis),dtype='str')
    subs = list(subs)
    print(len(subs),subs)

    for hemi in ['lh','rh']:
        for label in labels:
            try:
                saveProbLabelForAll(subs,hemi,label,outdir)
                saveProbLabelAverage(subs,hemi,label,outdir)
            except:
                pass

