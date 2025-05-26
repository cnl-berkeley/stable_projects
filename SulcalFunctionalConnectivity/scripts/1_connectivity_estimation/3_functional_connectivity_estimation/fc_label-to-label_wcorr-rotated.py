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


def procTask(postproc_dir,sub,ses,task,runs,subtask,labels,label_inds_all_lh,label_inds_all_rh,savedir,savename,highres_funcs_lh,highres_funcs_rh):
    n_labels=len(labels)

    # fix rh rot label indexing
    ts,affine = getRunData(highres_funcs_lh[0])
    for i in np.arange(len(label_inds_all_rh)):
        tmp = np.array(label_inds_all_rh[i]) - ts.shape[1]
        label_inds_all_rh[i] = list(tmp)

    y1_lh,w1,affine = procTaskRun(sub,ses,task,runs[0],subtask,highres_funcs_lh[0],label_inds_all_lh)
    y2_lh,w2,affine = procTaskRun(sub,ses,task,runs[1],subtask,highres_funcs_lh[1],label_inds_all_lh)
    y3_lh,w3,affine = procTaskRun(sub,ses,task,runs[2],subtask,highres_funcs_lh[2],label_inds_all_lh)

    y1_rh,w1,affine = procTaskRun(sub,ses,task,runs[0],subtask,highres_funcs_rh[0],label_inds_all_rh)
    y2_rh,w2,affine = procTaskRun(sub,ses,task,runs[1],subtask,highres_funcs_rh[1],label_inds_all_rh)
    y3_rh,w3,affine = procTaskRun(sub,ses,task,runs[2],subtask,highres_funcs_rh[2],label_inds_all_rh)

    y = []
    for l in arange(n_labels):
        y.append(concatenate((y1_lh[l],y2_lh[l],y3_lh[l]),axis=0))
    for l in arange(n_labels):
        y.append(concatenate((y1_rh[l],y2_rh[l],y3_rh[l]),axis=0))
    w = concatenate((w1,w2,w3),axis=0)

    stat_map = zeros((2*n_labels,2*n_labels))
    for i in arange(2*n_labels):
        for j in arange(2*n_labels):
            try:
                val = corr(y[i],y[j],w)
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
    print(ts.shape)

    y = []
    for label_inds in label_inds_all:
        if label_inds != []:
            x = mean(ts[:,label_inds], axis=1) # mean within label
        else:
            x = []
        y.append(x)

    if subtask[0] == 'hybrid':
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


def procSub(postproc_dir,fmriprep_derivatives_dir,sub,ses,task,subtasks,runs,label_names,rot_dir,tmpdir,outdir):
    savedir = '%s/sub-%s/ses-%s'%(outdir,sub,ses)
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    # collect lh functional data
    hemi = 'lh'
    highres_funcs_lh = []
    for run in runs:
        highres_func = '%s/sub-%s/ses-%s/hemi-%s/lev1/run-%s/res4d.func.gii'%(postproc_dir,sub,ses,hemi,run)
        highres_funcs_lh.append(highres_func)

    # collect rh functional data
    hemi = 'rh'
    highres_funcs_rh = []
    for run in runs:
        highres_func = '%s/sub-%s/ses-%s/hemi-%s/lev1/run-%s/res4d.func.gii'%(postproc_dir,sub,ses,hemi,run)
        highres_funcs_rh.append(highres_func)

    # build a rotated sulcal network and compute connectivity
    n_reps = 1000
    perm_inds = np.zeros((n_labels,n_reps))
    for l in np.arange(n_labels):
        inds = np.arange(n_reps)
        np.random.shuffle(inds)
        perm_inds[l,:] = inds

    for i in np.arange(n_reps):
        print('permutation %d'%i)
        hemi = 'lh'
        label_inds_all_lh = []
        for l in np.arange(n_labels):
            labelname = '%s.%s'%(hemi,labels[l])
            rot_file = '%s/sub-%s/rotated_labels/%s_rotations.csv'%(rot_dir,sub,labelname)
            label_inds_all_lh.append(getRotLabelInds(rot_file,perm_inds[l,i]))

        hemi = 'rh'
        label_inds_all_rh = []
        for l in np.arange(n_labels):
            labelname = '%s.%s'%(hemi,labels[l])
            rot_file = '%s/sub-%s/rotated_labels/%s_rotations.csv'%(rot_dir,sub,labelname)
            label_inds_all_rh.append(getRotLabelInds(rot_file,perm_inds[l,i]))

        if sub == 'n037t1':
            outdir = '%s/sub-%s/rotated_labels/labelfiles_iter%d'%(rot_dir,sub,i+1)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            saveLabelFiles(label_inds_all_lh,label_inds_all_rh,label_names,outdir) # testing only!!! subject hardcoded


def saveLabelFiles(label_inds_all_lh,label_inds_all_rh,label_names,save_dir):
    n_labels=len(labels)

    surf_file_lh = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7/sub-n037t1/surf/lh.inflated'
    surf_file_rh = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7/sub-n037t1/surf/rh.inflated'

    # fix rh rot label indexing
    for i in np.arange(len(label_inds_all_rh)):
        tmp = np.array(label_inds_all_rh[i]) - ts.shape[1]
        label_inds_all_rh[i] = list(tmp)

    for i in np.arange(len(label_names)):
        hemi = 'lh'
        tmp = nib.freesurfer.io.read_geometry(surf_file_lh)
        coords = tmp[0]
        vertices = tmp[1]

        out_file = '%s/lh.%s_ex_iteration.label'%(outdir,label_names[i])
        of = open(out_file,'w')
        of.write('# %s label for example iteration of null network\n'%label_names[i])
        n_vertices = len(label_inds_all_lh[i])
        of.write('%d\n'%n_vertices)
        for n in arange(n_vertices):
             vind = label_inds_all_lh[i][n]
             of.write('%d %.4f %.4f %.4f 0.00000\n'%(vind,coords[vind,0],coords[vind,1],coords[vind,2]))
        of.close()


        hemi = 'rh'
        tmp = nib.freesurfer.io.read_geometry(surf_file_rh)
        coords = tmp[0]
        vertices = tmp[1]

        out_file = '%s/rh.%s_ex_iteration.label'%(outdir,label_names[i])
        of = open(out_file,'w')
        of.write('# %s label for example iteration of null network\n'%label_names[i])
        n_vertices = len(label_inds_all_rh[i])
        of.write('%d\n'%n_vertices)
        for n in arange(n_vertices):
             vind = label_inds_all_rh[i][n]
             of.write('%d %.4f %.4f %.4f 0.00000\n'%(vind,coords[vind,0],coords[vind,1],coords[vind,2]))
        of.close()


def procLabels(postproc_dir,sub,ses,task,subtasks,labels,label_inds_all_lh,label_inds_all_rh,savedir,highres_funcs_lh,highres_funcs_rh,perm_ind):
    # loop wanted task conditions and hemispheres
    for subtask in subtasks:
        print('* subtask',subtask)
        savename = 'labels_%s_perm%d'%(subtask[0],perm_ind+1)
        procTask(postproc_dir,sub,ses,task,runs,subtask,labels,label_inds_all_lh,label_inds_all_rh,savedir,savename,highres_funcs_lh,highres_funcs_rh)


def getRotLabelInds(rot_file,rot_ind):
    line = np.loadtxt(rot_file,skiprows=int(rot_ind),max_rows=1,dtype='str')
    inds = str(line).strip('\n').split(',')
    inds = np.array(inds)
    inds = inds[inds != '']
    inds = np.array(inds).astype(int)
    inds = list(inds)

    return inds


if __name__ == '__main__':

    sub = sys.argv[1]
    ses = sub[-1]
    print(sub,ses)

    fmriprep_derivatives_dir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/bids/derivatives_v7'

    task = 'relmatch'
    tmpdir = '/home/weiner/shakki/NotBackedUp/bp-rot'

    postproc_dir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/postproc_fmri_smooth0mm_v7-bp/%s'%task
    outdir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2seed-bp/%s/rotated'%task

    lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p','prts','lfms','aalf']
    lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    labels = np.array(lpfc_labels + lpar_labels)
    n_labels = len(labels)

    runs = ['01','02','03']

    subtasks = [['hybrid',['hybrid']]]

    os.environ['SUBJECTS_DIR'] = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7'

    rot_dir = '/home/weiner/shakki/scripts-github/analysis/misc/generate_rotated_labels'

    procSub(postproc_dir,fmriprep_derivatives_dir,sub,ses,task,subtasks,runs,labels,rot_dir,tmpdir,outdir)

