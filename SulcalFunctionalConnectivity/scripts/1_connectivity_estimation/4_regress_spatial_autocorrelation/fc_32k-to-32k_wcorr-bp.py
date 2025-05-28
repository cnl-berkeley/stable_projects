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
    giidata = nib.load(labelfile)
    label_vec = squeeze(asarray([x.data for x in giidata.darrays]))
    label_inds = where(label_vec > 0)[0]
    label_inds = label_inds.astype(int32)
#    im_label = nib.freesurfer.read_label(labelfile)
    return list(label_inds)


def procTask(postproc_dir,sub,ses,hemi,task,runs,subtask,label_inds,outdir,savename,lowres_funcs):

    y1,w1,affine = procTaskRun(sub,ses,hemi,task,runs[0],subtask,label_inds,lowres_funcs[0])
    y2,w2,affine = procTaskRun(sub,ses,hemi,task,runs[1],subtask,label_inds,lowres_funcs[1])
    y3,w3,affine = procTaskRun(sub,ses,hemi,task,runs[2],subtask,label_inds,lowres_funcs[2])

    y = concatenate((y1,y2,y3),axis=0)
    w = concatenate((w1,w2,w3),axis=0)

    n_points=shape(y)[1]

    y = y[w > 0,:]
    stat_map = corrcoef(y,rowvar=False)
    stat_map_z = arctanh(stat_map)

    print(stat_map_z)
    stat_map_z = stat_map.astype('float32')

    stat_map_z = around(stat_map_z,decimals=3) # save three decimals

    savez_compressed('%s/fc-%s-%s.npz'%(outdir,savename,hemi),stat_map_z)


def procTaskRun(sub,ses,hemi,task,run,subtask,label_inds,lowres_func):
    y,affine = getRunData(lowres_func)
    y = y[:,label_inds]

    if subtask[0] == 'hybrid':
        w = ones(len(y))
    else:
        w = getRunTimingsWeighted(postproc_dir,sub,ses,run,hemi,task,subtask,len(y))
    print('w',w)

    w = zeroOutliers(w,postproc_dir,sub,ses,run)

    print('w_final',len(w),w)

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


def getRunTimingsWeighted(postproc_dir,sub,ses,run,hemi,task,subtask,n_scans):
    dropvols=3
    tr=2
    bids_dir='/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/bids'
    print('%s/*task-%s_run-%s_events.tsv'%(bids_dir,task,int(run)))
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

    return w


def getRunData(highres_func):
    im_file=nib.load(highres_func)
    y=squeeze(asarray([x.data for x in im_file.darrays]))
    # affine info is lost at postprocessing with workbench,
    # https://www.mail-archive.com/hcp-users@humanconnectome.org/msg06178.html
    affine = eye(4)
    return y, affine


def createSurfaces(fmriprep_derivatives_dir,sub,hemi,outdir):

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

    highres_inflated = '%s/sub-%s/ses-%s/anat/sub-%s_ses-%s_hemi-%s_inflated.surf.gii'%(fmriprep_derivatives_dir,sub,ses,sub,ses,hemi[0].upper())
    lowres_inflated = '%s/%s-32k.surf.gii'%(outdir,os.path.basename(highres_inflated).split('.')[0])
    str = 'wb_command -surface-resample %s %s %s BARYCENTRIC %s'%(highres_inflated,highres_sphere,lowres_sphere,lowres_inflated)
    cmd=wb.WBCommand(str)
    cmd.run()

    return lowres_sphere,highres_sphere,lowres_surface,lowres_inflated


def saveSurfacesFS(lowres_sphere,lowres_mid,lowres_inflated,outdir):
    import nipype.interfaces.freesurfer as fs
    mris = fs.MRIsConvert()
    for s in [lowres_sphere,lowres_mid,lowres_inflated]:    
        mris.inputs.in_file = s
        mris.inputs.out_file = '%s/%s'%(outdir,'.'.join(os.path.basename(s).split('.')[0:-2]))
        mris.run()


def downSampleFunc(hemi,highres_func,lowres_sphere,highres_sphere,lowres_surface,lowres_func):
    hemiletter = hemi[0].upper()

    highres_surface = '%s/sub-%s/ses-%s/anat/sub-%s_ses-%s_hemi-%s_midthickness.surf.gii'%(fmriprep_derivatives_dir,sub,ses,sub,ses,hemiletter)
    str = 'wb_command -metric-resample %s %s %s ADAP_BARY_AREA -area-surfs %s %s %s'%(highres_func,highres_sphere,lowres_sphere,highres_surface,lowres_surface,lowres_func)
    cmd=wb.WBCommand(str)
    cmd.run()

    return lowres_func


def downSampleLabel(hemi,labelname,lowres_sphere,highres_sphere,lowres_surface,outdir,fmriprep_derivatives_dir):
    # resample label to 32k space

    highres_label_bin = '%s/sub-%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,labelname)
 
    # convert to gii
    highres_surface_bin = '%s/sub-%s/surf/%s.white'%(os.environ['SUBJECTS_DIR'],sub,hemi)
    highres_label_bin = '%s/sub-%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,labelname)
    highres_label = '%s/%s.%s.label.gii'%(outdir,labelname,hemi[0].upper())
    str = 'mris_convert --label %s %s %s %s'%(highres_label_bin,labelname,highres_surface_bin,highres_label)
    os.system(str)

    # resample
    highres_surface = '%s/sub-%s/ses-%s/anat/sub-%s_ses-%s_hemi-%s_midthickness.surf.gii'%(fmriprep_derivatives_dir,sub,ses,sub,ses,hemi[0].upper())
    lowres_label = '%s/%s-32k.%s.label.gii'%(outdir,labelname,hemi[0].upper())
    import nipype.interfaces.workbench.base as wb
    str = 'wb_command -label-resample %s %s %s ADAP_BARY_AREA -area-surfs %s %s %s'%(highres_label,highres_sphere,lowres_sphere,highres_surface,lowres_surface,lowres_label)
    cmd = wb.WBCommand(str)
    cmd.run()

    return lowres_label


def procSub(postproc_dir,fmriprep_derivatives_dir,sub,ses,task,subtasks,runs,labels,tmpdir,outdir):
    for hemi in ['lh','rh']:
        savedir = '%s/sub-%s/ses-%s/hemi-%s'%(outdir,sub,ses,hemi)
        if not os.path.exists(savedir):
            os.makedirs(savedir)

        # crate lowres surfaces to downsample functional data
        savedir_surf = '%s/sub-%s/surf_32k/gifti'%(outdir,sub)
        if not os.path.exists(savedir_surf):
            os.makedirs(savedir_surf)
        lowres_sphere,highres_sphere,lowres_surface,lowres_inflated = createSurfaces(fmriprep_derivatives_dir,sub,hemi,savedir_surf)
        # save lowres surfaces in FreeSurfer format for later plotting
        savedir_surf = '%s/sub-%s/surf_32k/fs'%(outdir,sub)
        if not os.path.exists(savedir_surf):
            os.makedirs(savedir_surf)
        saveSurfacesFS(lowres_sphere,lowres_surface,lowres_inflated,savedir_surf)

        # downsample functional data
        tmp_savedir = '%s/sub-%s/ses-%s/hemi-%s'%(tmpdir,sub,ses,hemi)
        if not os.path.exists(tmp_savedir):
            os.makedirs(tmp_savedir)

        highres_funcs = [] # to extract label timecourses
        lowres_funcs = [] # brain parcellation
        for run in runs:
            highres_func = '%s/sub-%s/ses-%s/hemi-%s/lev1/run-%s/res4d.func.gii'%(postproc_dir,sub,ses,hemi,run)
            highres_funcs.append(highres_func)
            lowres_savepath = '%s/sub-%s_ses-%s_hemi-%s_run-%s_res4d-32k.func.gii'%(tmp_savedir,sub,ses,hemi,run)
            lowres_func = downSampleFunc(hemi,highres_func,lowres_sphere,highres_sphere,lowres_surface,lowres_savepath)
            lowres_funcs.append(lowres_func)

        for label in labels:
            # get seed label indices (lowres space)
            lowres_label = downSampleLabel(hemi,label,lowres_sphere,highres_sphere,lowres_surface,tmp_savedir,fmriprep_derivatives_dir)
            label_inds = getLabelInds(lowres_label)
            procLabel(postproc_dir,sub,ses,hemi,task,subtasks,label,label_inds,savedir,lowres_funcs)


def procLabel(postproc_dir,sub,ses,hemi,task,subtasks,label,label_inds,savedir,lowres_funcs):
    # loop wanted task conditions and hemispheres
    for subtask in subtasks:
        print('* subtask',subtask)
        savename = '%s_%s_%s'%(hemi,label,subtask[0])
        procTask(postproc_dir,sub,ses,hemi,task,runs,subtask,label_inds,savedir,savename,lowres_funcs)


if __name__ == '__main__':
    fmriprep_derivatives_dir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/bids/derivatives_v7'

    task = 'relmatch'
    postproc_dir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/postproc_fmri_smooth0mm_v7-bp/%s'%task

    labels = ['cortex'] #['cortex']

    subtasks = [['hybrid',['hybrid']]]

    os.environ['SUBJECTS_DIR'] = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7'

    sub = sys.argv[1]
    ses = sub[-1]
    runs = ['01','02','03']
    print(sub,ses,runs)

    tmpdir='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation_tmp'
    outdir='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/cortex/fc_32k232k_0mm-bp-%s'%(sub,task)
    procSub(postproc_dir,fmriprep_derivatives_dir,sub,ses,task,subtasks,runs,labels,tmpdir,outdir)

