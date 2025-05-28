import nibabel as nib
from numpy import *
import os,sys
import pandas as pd
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
        im_label = nib.freesurfer.read_label(labelfile)
    except:
        im_label = []

    return list(im_label)


def getAnnotInds(annot_file):
    labels,ctab,names = nib.freesurfer.io.read_annot(annot_file)
    uni_vals = unique(labels) 
    annot_inds = []
    names_str = []
    for i in arange(len(uni_vals)):
        name = names[i] #.decode('utf-8')
        if name != 'unknown':
            annot_inds.append(list(where(labels == uni_vals[i])[0]))

            if name not in names_str:
                names_str.append(name)

    return names_str,annot_inds


def procTask(postproc_dir,sub,ses,hemi,task,runs,subtask,labels,label_inds_all,annots,annot_inds,outdir,savename,highres_funcs):
    # labels
    n_labels=len(label_inds_all)
    y1,w1,affine = procTaskRun(sub,ses,hemi,task,runs[0],subtask,highres_funcs[0],label_inds_all)
    y2,w2,affine = procTaskRun(sub,ses,hemi,task,runs[1],subtask,highres_funcs[1],label_inds_all)
    y3,w3,affine = procTaskRun(sub,ses,hemi,task,runs[2],subtask,highres_funcs[2],label_inds_all)
    y = []
    for l in arange(n_labels):
        y.append(concatenate((y1[l],y2[l],y3[l]),axis=0))
    w = concatenate((w1,w2,w3),axis=0)

    # annots
    n_annots = len(annot_inds)
    y1b,w1,affine = procTaskRun(sub,ses,hemi,task,runs[0],subtask,highres_funcs[0],annot_inds)
    y2b,w2,affine = procTaskRun(sub,ses,hemi,task,runs[1],subtask,highres_funcs[1],annot_inds)
    y3b,w3,affine = procTaskRun(sub,ses,hemi,task,runs[2],subtask,highres_funcs[2],annot_inds)
    yb = []
    for l in arange(n_annots):
        yb.append(concatenate((y1b[l],y2b[l],y3b[l]),axis=0))

    stat_map = zeros((n_labels,n_annots))
    for i in arange(n_labels):
        seed_ts = y[i]
        for j in arange(n_annots):
            #print(i,j)
            annot_ts = yb[j]
            val = corr(seed_ts,annot_ts,w)
            stat_map[i,j] = arctanh(val)
    stat_map = stat_map.astype('float32')
    outpath = '%s/seed-fc-annot_%s.csv'%(outdir,savename)
    df = pd.DataFrame(stat_map,columns=annots)
    df.to_csv(outpath,index=False)


def procTaskRun(sub,ses,hemi,task,run,subtask,highres_func,label_inds_all):
    y = []
    for label_inds in label_inds_all:
        ts,affine = getRunData(highres_func)
        x = mean(ts[:,label_inds], axis=1) # mean within label
        y.append(x)

    if subtask[0] == 'hybrid':
        w = ones(len(x))
    else:
        w = getRunTimingsWeighted(postproc_dir,sub,ses,run,hemi,task,subtask,len(x))

    return y,w,affine


def getRunTimingsWeighted(postproc_dir,sub,ses,run,hemi,task,subtask,n_scans):
    dropvols=3
    tr=2
    bids_dir='/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/bids'
    events_file = glob('%s/*task-%s_run-%s_events.tsv'%(bids_dir,task,int(run)))[0]
    print(events_file)
    event_df = pd.read_csv(events_file,delim_whitespace=True)
    condition_names = list(unique(event_df.trial_type))
    frame_times = arange(0,tr*n_scans,tr)
    w = [0]*n_scans
    for cond in subtask[1]:
        onsets = list(event_df[event_df.trial_type == cond].onset - 1.0*tr*dropvols)
        durations = list(event_df[event_df.trial_type == cond].duration.astype(float))
        if cond == 'rest':
            onsets = list(event_df[event_df.trial_type == 'instruction'].onset - 1.0*tr*dropvols)
            onsets = [x-14 for x in onsets]
            durations = [14.0]*len(onsets)
        amplitudes = [1]*len(onsets)

        w_cond,w_name = compute_regressor((onsets,durations,amplitudes),'spm',frame_times)
        w_cond[w_cond < 0] = 0 # only accept positive weights
        w_cond = list(squeeze(w_cond))
        w = [max(l1,l2) for l1,l2 in zip(w,w_cond)]

    #  exclude outliers
    outliers_ind_file = '%s/sub-%s/ses-%s/hemi-%s/lev1/run-%s/event-inds_outliers.csv'%(postproc_dir,sub,ses,hemi,run) # zero-based indices
    try:
        outlier_inds = list(squeeze(array(pd.read_csv(outliers_ind_file,header=None,sep='\n'))))
        w[outlier_inds] = 0
    except:
       # assume no outliers
       pass

    return w


def getRunData(highres_func):
    im_file=nib.load(highres_func)
    y=squeeze(asarray([x.data for x in im_file.darrays]))
    # affine info is lost at postprocessing with workbench,
    # https://www.mail-archive.com/hcp-users@humanconnectome.org/msg06178.html
    affine = ones(4)
    return y, affine


def procSub(postproc_dir,fmriprep_derivatives_dir,sub,ses,task,subtasks,runs,labels,tmpdir,outdir):
    for hemi in ['lh','rh']:
        savedir = '%s/sub-%s/ses-%s/hemi-%s'%(outdir,sub,ses,hemi)
        if not os.path.exists(savedir):
            os.makedirs(savedir)

        highres_funcs = []
        for run in runs:
            highres_func = '%s/sub-%s/ses-%s/hemi-%s/lev1/run-%s/res4d.func.gii'%(postproc_dir,sub,ses,hemi,run)
            highres_funcs.append(highres_func)

        label_inds_all = []
        for label in labels:
            if type(label) is list: # pick the first label in list that exits
                label_found = False
                for lab in label:
                    labelfile = '%s/sub-%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,lab)
                    if os.path.exists(labelfile) and not label_found:
                        label_inds = getLabelInds(labelfile)
                        label_inds_all.append(label_inds)
                        label_found = True
                if label_found == False:
                    label_inds_all.append([])
            else: # look for specific label
                labelfile = '%s/sub-%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,label)
                label_inds = getLabelInds(labelfile)
                label_inds_all.append(label_inds)

        procLabels(postproc_dir,sub,ses,hemi,task,subtasks,labels,label_inds_all,savedir,highres_funcs)


def procLabels(postproc_dir,sub,ses,hemi,task,subtasks,labels,label_inds_all,savedir,highres_funcs):
    annot_file='/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/seed-fc/mibd/annot-relmatch/sub-%s_ses-%s_combined_clusters_thresh0.61-%s-nw-full.annot'%(sub,sub[-1],hemi)
    annots,annot_inds = getAnnotInds(annot_file)

    # loop wanted task conditions and hemispheres
    for subtask in subtasks:
        print('* subtask',subtask)
        savename = '%s_nw_%s'%(hemi,subtask[0])
        procTask(postproc_dir,sub,ses,hemi,task,runs,subtask,labels,label_inds_all,annots,annot_inds,savedir,savename,highres_funcs)


if __name__ == '__main__':
    fmriprep_derivatives_dir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/bids/derivatives_v7'

    task = 'relmatch'
    postproc_dir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/postproc_fmri_smooth0mm_v7-bp/%s'%task

    tmpdir='/home/weiner/shakki/NotBackedUp'
    outdir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2annot-bp/%s'%task

    lpfc_labels = ['ifs',['painfs_any','painfs_combined'],'pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
    lpar_labels = [['slos1','slocs-v','SLOS'],'sB','pips','mTOS',['iTOS','ITOS','lTOS'],'IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    mpar_labels = []

    labels = lpfc_labels + lpar_labels + mpar_labels

    runs = ['01','02','03']

    subtasks = [['hybrid',['hybrid']]]

    os.environ['SUBJECTS_DIR'] = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7'

    sub = sys.argv[1]
    ses = sub[-1]
    print(sub,ses,runs)

    procSub(postproc_dir,fmriprep_derivatives_dir,sub,ses,task,subtasks,runs,labels,tmpdir,outdir)

