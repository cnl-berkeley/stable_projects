import os,sys
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp

def getLabelInd(labels,label):
    return np.argwhere(labels == label)[0][0]


def getOrigData(datadir,sub,hemi,labels,label):
    fc_file = '%s/sub-%s/ses-%s/hemi-%s/seed-fc-annot_%s_nw_hybrid.csv'%(datadir,sub,sub[-1],hemi,hemi)
    print(fc_file)
    data = pd.read_csv(fc_file)
    column_names = data.columns
    data = np.array(data)[getLabelInd(labels,label),:]
    data = data[1::] # exclude 'unknown' label
    print(np.shape(data))
    column_names = np.array(column_names[1::])
    for c in np.arange(len(column_names)):
       column_names[c] = column_names[c][2:-1]

    return data,column_names


def getRotData(datadir,sub,hemi,label):
    fc_file = '%s/sub-%s/ses-%s/hemi-%s/seed-fc-annot_rotated-%s_%s.%s_nw_hybrid.csv'%(datadir,sub,sub[-1],hemi,hemi,hemi,label)
    print(fc_file)
    data = pd.read_csv(fc_file)
    column_names = data.columns
    data = np.array(data)
    data = data[:,1::] # exclude 'unknown' label
    print(np.shape(data))
    column_names =  np.array(column_names[1::])
    for c in np.arange(len(column_names)):
       column_names[c] = column_names[c][2:-1]

    return data,column_names


def doOneSulcus(datadir,subs,labels,label,hemi,outdir):
    # collect data
    n_reps = 1000
    n_nw = 14
    all_orig = np.zeros((n_subs,n_nw))
    all_rot = np.zeros((n_subs,n_nw,n_reps))
    for s in np.arange(n_subs):
        orig,nwnames = getOrigData(datadir,subs[s],hemi,labels,label)
        all_orig[s,:] = orig

        rot,nwnames = getRotData(datadir,subs[s],hemi,label)
        all_rot[s,:,:] = rot.T

    # bootstrap without replacement
    # at each iteration use paired t-test
    perm_inds = np.zeros((n_subs,n_reps))
    for s in np.arange(n_subs):
        inds = np.arange(n_reps)
        np.random.shuffle(inds)
        perm_inds[s,:] = inds

    # is sulcal correlation significant different from random patches?
    collected_results = np.zeros((n_nw,4))*np.nan # orig t, rot t, p_right, p_left
    for nw in np.arange(n_nw):
        res = ttest_1samp(all_orig[:,nw],0)
        orig_t = res.statistic
        
        rot_ts = np.zeros((n_reps,1))*np.nan
        for i in np.arange(n_reps):
            rot_data = np.zeros((n_subs,1))*np.nan
            for s in np.arange(n_subs):
                rot_data[s] = all_rot[s,nw,perm_inds[s,:] == i]
            rot_res = ttest_1samp(rot_data,0)
            rot_ts[i] = rot_res.statistic[0]

        collected_results[nw,0] = orig_t
        collected_results[nw,1] = np.mean(rot_ts)
        boot_p = len(np.argwhere(rot_ts > orig_t))/n_reps # right tail (original more correlated than rotated)
        collected_results[nw,2] = boot_p

        boot_p = len(np.argwhere(rot_ts < orig_t))/n_reps # left tail (original less correlated than rotated)
        collected_results[nw,3] = boot_p

    np.savez('%s/collected_results_%s.%s.npz'%(outdir,hemi,label),collected_results,nwnames)


if __name__ == '__main__':

    outdir = '%s/bootstrap'%os.getcwd()
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p','prts','lfms','aalf']
    lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    labels = np.array(lpfc_labels + lpar_labels)

    datadir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2annot-bp/relmatch'

    namelistfile = '/home/weiner/shakki/scripts-new/preprocessing/regress_motion/-bp/orig/fc-sareg-fdreg/namelist.csv'
    subs = np.loadtxt(namelistfile,dtype='str')
    n_subs = len(subs)
    print('n_subs',n_subs,subs)

    for hemi in ['lh','rh']:
        for label in labels:
            try:
                doOneSulcus(datadir,subs,labels,label,hemi,outdir)
            except:
                pass

