import nibabel as nib
from numpy import *
import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import nipype.interfaces.workbench.base as wb
from nilearn.glm.first_level import compute_regressor
from glob import glob

def m(x, w):
    """Weighted Mean"""
    return sum(x * w) / sum(w)


def cov(x, y, w):
    """Weighted Covariance"""
    return sum(w * (x - m(x, w)) * (y - m(y, w))) / sum(w)


def corr(x, y, w):
    """Weighted Correlation"""
    # https://stackoverflow.com/questions/64468816/weighted-correlation-in-python
    return cov(x, y, w) / sqrt(cov(x, x, w) * cov(y, y, w))


def getLabelInds(labelfile):
    try:
        im_label = list(nib.freesurfer.read_label(labelfile))
    except:
        im_label = []

    return im_label


def procTask(postproc_dir,sub,ses,task,runs,subtask,labels,label_inds_all_lh,label_inds_all_rh,savedir,savename,highres_funcs_lh,highres_funcs_rh):
    n_labels=len(labels)

    y1_lh,w1,affine = procTaskRun(sub,ses,task,runs[0],subtask,highres_funcs_lh[0],label_inds_all_lh)
    y2_lh,w2,affine = procTaskRun(sub,ses,task,runs[1],subtask,highres_funcs_lh[1],label_inds_all_lh)
    y3_lh,w3,affine = procTaskRun(sub,ses,task,runs[2],subtask,highres_funcs_lh[2],label_inds_all_lh)

    y1_rh,w1,affine = procTaskRun(sub,ses,task,runs[0],subtask,highres_funcs_rh[0],label_inds_all_rh)
    y2_rh,w2,affine = procTaskRun(sub,ses,task,runs[1],subtask,highres_funcs_rh[1],label_inds_all_rh)
    y3_rh,w3,affine = procTaskRun(sub,ses,task,runs[2],subtask,highres_funcs_rh[2],label_inds_all_rh)

    print(len(highres_funcs_lh),len(highres_funcs_rh))
    print(len(label_inds_all_lh),len(label_inds_all_rh),n_labels,labels)
    y = []
    for l in arange(n_labels):
        y.append(concatenate((y1_lh[l],y2_lh[l],y3_lh[l]),axis=0))
    for l in arange(n_labels):
        y.append(concatenate((y1_rh[l],y2_rh[l],y3_rh[l]),axis=0))
    w = concatenate((w1,w2,w3),axis=0)

    print('y',len(y),y)
    print('w',len(w),w)
    print('n_labels',n_labels,labels)

    stat_map = zeros((2*n_labels,2*n_labels))
    for i in arange(2*n_labels):
        for j in arange(2*n_labels):
            try:
                val = corr(y[i],y[j],w)
                #val = corrcoef(y[i],y[j])[0,1]
                print('corr',val)
                stat_map[i,j] = arctanh(val)
            except:
                stat_map[i,j] = NaN

    print('stat_map',stat_map)

    stat_map = stat_map.astype('float32')

    column_names = []
    for hemi in ['lh','rh']:
        for label in labels:
            column_names.append('%s.%s'%(hemi,label))

    outpath = '%s/seed-fc_%s.csv'%(savedir,savename)
    print('saving outpath',outpath)
    df = pd.DataFrame(stat_map,columns=column_names)
    df.to_csv(outpath,index=False)


def procTaskRun(sub,ses,task,run,subtask,highres_func,label_inds_all):
    ts,affine = getRunData(highres_func)

    y = []
    for label_inds in label_inds_all:
        print('label_inds',len(label_inds))
        if label_inds != []:
            x = mean(ts[:,label_inds], axis=1) # mean within label
        else:
            x = []
        print('len of x',len(x))
        y.append(x)

    if subtask[0] == 'hybrid':
        print('len w', len(ts[:,1]))
        w = ones(len(ts[:,1]))
    else:
        w = getRunTimingsWeighted(postproc_dir,sub,ses,run,task,subtask,len(ts[:,1]))

    w = zeroOutliers(w,postproc_dir,sub,ses,run)

    return y,w,affine


def zeroOutliers(w,postproc_dir,sub,ses,run):
    #  exclude outliers
    outliers_ind_file = '%s/sub-%s/ses-%s/hemi-%s/lev1/run-%s/event-inds_outliers.csv'%(postproc_dir,sub,ses,'lh',run) # zero-based indices
    try:
        outlier_inds = list(squeeze(array(pd.read_csv(outliers_ind_file,header=None,sep='\n'))))
        w[outlier_inds] = 0
    except:
       # assume no outliers
       pass

    return w


def getRunTimingsWeighted(postproc_dir,sub,ses,run,task,subtask,n_scans):
    dropvols=3
    tr=2
    bids_dir='/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/bids'
    print('%s/*task-%s_run-%s_events.tsv'%(bids_dir,task,int(run)))
    events_file = glob('%s/*task-%s_run-%s_events.tsv'%(bids_dir,task,int(run)))[0]
    print(events_file)
    event_df = pd.read_csv(events_file,delim_whitespace=True)
    frame_times = arange(0,tr*n_scans,tr)
    w = [0]*n_scans
    for cond in subtask[1]:
        print('cond',cond)
        onsets = list(event_df[event_df.trial_type == cond].onset - 1.0*tr*dropvols)
        durations = list(event_df[event_df.trial_type == cond].duration.astype(float))
        if cond == 'rest':
            onsets = list(event_df[event_df.trial_type == 'instruction'].onset - 1.0*tr*dropvols)
            onsets = [x-14 for x in onsets]
            durations = [14.0]*len(onsets)
        amplitudes = [1]*len(onsets)

        print('onsets',onsets)
        print('durations',durations)
        try:
            w_cond,w_name = compute_regressor((onsets,durations,amplitudes),'spm',frame_times)
            w_cond[w_cond < 0] = 0 # only accept positive weights
            w_cond = list(squeeze(w_cond))
            w = [max(l1,l2) for l1,l2 in zip(w,w_cond)]
        except:
            w = []
    return w


def getRunData(highres_func):
    im_file=nib.load(highres_func)
    y=squeeze(asarray([x.data for x in im_file.darrays]))
    # affine info is lost at postprocessing with workbench,
    # https://www.mail-archive.com/hcp-users@humanconnectome.org/msg06178.html
    affine = ones(4)
    return y, affine


def createSurfaces(fmriprep_derivatives_dir,sub,hemi,outdir):
    # upsample functional
    current_sphere_bin = '%s/sub-%s/surf/%s.sphere.reg'%(os.environ['SUBJECTS_DIR'],sub,hemi)
    highres_sphere = '%s/%s.sphere.reg.surf.gii'%(outdir,hemi)
    import nipype.interfaces.freesurfer as fs
    mris = fs.MRIsConvert()
    mris.inputs.in_file = current_sphere_bin
    mris.inputs.out_file = highres_sphere
    mris.run()

    # create 32k mesh (2mm spacing)
    # example in https://www.humanconnectome.org/software/workbench-command/-surface-create-sphere
    lowres_sphere = '%s/sphere.32k.%s.surf.gii'%(outdir,hemi[0].upper())
    str = 'wb_command -surface-create-sphere 32000 %s'%lowres_sphere
    cmd=wb.WBCommand(str)
    cmd.run()

    if hemi == 'rh':
        str = 'wb_command -set-structure %s CORTEX_RIGHT'%lowres_sphere
        cmd=wb.WBCommand(str)
        cmd.run()
    else:
        str = 'wb_command -surface-flip-lr %s %s'%(lowres_sphere,lowres_sphere)
        cmd=wb.WBCommand(str)
        cmd.run()

        str = 'wb_command -set-structure %s CORTEX_LEFT'%lowres_sphere
        cmd=wb.WBCommand(str)
        cmd.run()

    highres_surface = '%s/sub-%s/ses-%s/anat/sub-%s_ses-%s_hemi-%s_midthickness.surf.gii'%(fmriprep_derivatives_dir,sub,ses,sub,ses,hemi[0].upper())
    lowres_surface = '%s/%s-32k.surf.gii'%(outdir,os.path.basename(highres_surface).split('.')[0])
    str = 'wb_command -surface-resample %s %s %s BARYCENTRIC %s'%(highres_surface,highres_sphere,lowres_sphere,lowres_surface)
    cmd=wb.WBCommand(str)
    cmd.run()

    return lowres_sphere,highres_sphere,lowres_surface


def findLabelFile(label,sub,hemi,prob_tag):
    if prob_tag == 0:
        labeldir = '%s/sub-%s/label'%(os.environ['SUBJECTS_DIR'],sub)
    else:
        labeldir = '%s/sub-%s/label/sulcal_prob_mpm'%(os.environ['SUBJECTS_DIR'],sub)

    if type(label) is list: # pick the first label in list that exits
        for lab in label:
            labelfile = '%s/%s.%s.label'%(labeldir,hemi,lab)
            if os.path.exists(labelfile):
                return labelfile
            else:
                labelfile = '%s/%s.%s.label'%(labeldir,hemi,lab)
                if os.path.exists(labelfile):
                    return labelfile
 
    else: # look for specific label
        labelfile = '%s/%s.%s.label'%(labeldir,hemi,label)
        if os.path.exists(labelfile):
            return labelfile
            labelfile = '%s/%s.%s.label'%(labeldir,hemi,label)
            if os.path.exists(labelfile):
                return labelfile
    print(labelfile)

    return ''


def procSub(postproc_dir,fmriprep_derivatives_dir,sub,ses,task,subtasks,runs,labels,tmpdir,outdir,prob_tag):
    savedir = '%s/sub-%s/ses-%s'%(outdir,sub,ses)
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    # collect lh functional data
    hemi = 'lh'
    highres_funcs_lh = []
    for run in runs:
        highres_func = '%s/sub-%s/ses-%s/hemi-%s/lev1/run-%s/res4d.func.gii'%(postproc_dir,sub,ses,hemi,run)
        highres_funcs_lh.append(highres_func)

    label_inds_all_lh = []
    label_names = []
    for label in labels:
        # get seed label indices (native space)
        labelfile = findLabelFile(label,sub,hemi,prob_tag)
        label_inds = getLabelInds(labelfile)
        label_inds_all_lh.append(label_inds)

        # store label names
        if type(label) is list:
            label_names.append(label[0])
        else:
            label_names.append(label)

    # collect rh functional data
    hemi = 'rh'
    highres_funcs_rh = []
    for run in runs:
        highres_func = '%s/sub-%s/ses-%s/hemi-%s/lev1/run-%s/res4d.func.gii'%(postproc_dir,sub,ses,hemi,run)
        highres_funcs_rh.append(highres_func)

    label_inds_all_rh = []
    for label in labels:
        # get seed label indices (native space)
        labelfile = findLabelFile(label,sub,hemi,prob_tag)
        label_inds = getLabelInds(labelfile)
        label_inds_all_rh.append(label_inds)

    procLabels(postproc_dir,sub,ses,task,subtasks,label_names,label_inds_all_lh,label_inds_all_rh,savedir,highres_funcs_lh,highres_funcs_rh)


def procLabels(postproc_dir,sub,ses,task,subtasks,labels,label_inds_all_lh,label_inds_all_rh,savedir,highres_funcs_lh,highres_funcs_rh):
    # loop wanted task conditions and hemispheres
    for subtask in subtasks:
        print('* subtask',subtask)
        savename = 'labels_%s'%subtask[0]
        procTask(postproc_dir,sub,ses,task,runs,subtask,labels,label_inds_all_lh,label_inds_all_rh,savedir,savename,highres_funcs_lh,highres_funcs_rh)


if __name__ == '__main__':

    sub = sys.argv[1]
    ses = sub[-1]
    print(sub,ses)

    fmriprep_derivatives_dir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/bids/derivatives_v7'

    task = 'relmatch'
    tmpdir = '/home/weiner/shakki/NotBackedUp/bp-prob'

    postproc_dir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/postproc_fmri_smooth0mm_v7-bp/%s'%task
    outdir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2seed-bp/%s'%task

    lpfc_labels = ['ifs','painfs_a','painfs_p',['painfs_any','painfs_combined'],'pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
    lpar_labels = [['STS','sts'],['SPS','sps'],['SmgS','smgS','Smgs'],['slos1','slocs-v','SLOS'],'sB','pips','mTOS',['iTOS','ITOS','lTOS'],'IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    mpar_labels = ['1','2','3','MCGS','POS','prculs','prcus1','prcus2','prcus3','sbps','sps']
    labels = lpfc_labels + lpar_labels + mpar_labels

    runs = ['01','02','03']

    subtasks = [['hybrid',['hybrid']]]

    os.environ['SUBJECTS_DIR'] = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7'

    procSub(postproc_dir,fmriprep_derivatives_dir,sub,ses,task,subtasks,runs,labels,tmpdir,outdir,0)


    ### repeat sulcal connectivity estimation using probabilistic labels

    lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p','prts','lfms','aalf']
    lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    labels = lpfc_labels + lpar_labels

    outdir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2seed-bp-prob_mpm/%s'%task

    procSub(postproc_dir,fmriprep_derivatives_dir,sub,ses,task,subtasks,runs,labels,tmpdir,outdir,1)

