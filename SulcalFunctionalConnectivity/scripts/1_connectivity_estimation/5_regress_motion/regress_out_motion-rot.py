import os,sys
import re
from numpy import *
import pandas as pd
from scipy import stats
import numpy as np

def getData(df):
    # the connectivity matrix is required to have non-nan value for all connections
    data_all = []
    sub_ids = []
    sub_fds = []
    for subindex, row in df.iterrows():
        sub = row.both_id
        ses = sub[-1]

        dfile = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2seed-bp/relmatch/rotated/sub-n002t2/ses-2/seed-fc_labels_hybrid_perm1000.csv'

        dfile = '/home/weiner/shakki/NotBackedUp/spatial_autocorrelation%s/%s/labelnw/fc-labels-both%s.txt'%(pipe,sub,tag)
        print(dfile)
        if os.path.exists(dfile):
            df_sub = pd.read_csv(dfile,sep='\s+',header=None)
            df_sub = df_sub.round(decimals=3)
            if not any(isnan(df_sub)):
                df_sub[isinf(df_sub)] = 0
                data_all.append(df_sub)
                sub_ids.append(sub)
                sub_fds.append(row['mFD'])

    return data_all,sub_fds,sub_ids
    

def getResiduals(data,correlate):
    print(data)
    n_subs = len(data)
    n_roi = shape(data[0])[1]
    x = correlate

    res_arr = zeros((n_subs,n_roi,n_roi))
    for i in arange(n_roi):
        for j in arange(n_roi):
            y = []
            for s in arange(n_subs):
                y.append(array(data[s])[i,j]) 
            
            ret = stats.linregress(x,y)
            for s in arange(n_subs):
               res_arr[s,i,j] = y[s] - ret[0]*x[s] # leave intercept

    return res_arr


def saveCohortData(data_all,sub_ids,savedir):
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    # save residualized data
    all_vals = zeros((len(sub_ids),len(data_all[0]),len(data_all[0])))
    for i in arange(len(sub_ids)):
        data = array(data_all[i])
        all_vals[i,:,:] = array(data)

    savez('%s/allsub-i%d-fd.npz'%(savedir,i),data=all_vals)


if __name__ == '__main__':

    sub_ids = np.loadtxt('/home/weiner/shakki/scripts-new/analysis/output/setup/-bp/orig/setup-sareg-fdreg/nonorm-both/regs/namelist.csv',dtype='str')
    fd_vals = np.loadtxt('/home/weiner/shakki/scripts-new/analysis/output/setup/-bp/orig/setup-sareg-fdreg/nonorm-both/regs/mFD.csv')
    datadir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2seed-bp/relmatch/rotated'

    print(len(fd_vals))

    n_subs = len(sub_ids)
    print(n_subs)
    n_reps = 2

    for i in np.arange(n_reps):
        data_all = []
        for s in np.arange(n_subs):
            dfile = '%s/sub-%s/ses-%s/seed-fc_labels_hybrid_perm%d.csv'%(datadir,sub_ids[s],sub_ids[s][-1],i+1)
            df_sub = pd.read_csv(dfile,sep='\s+',header=None)
            df_sub = df_sub.round(decimals=3)
            data_all.append(df_sub)
        
        # regress out effects of mfd
        res_all = getResiduals(data,fd_vals)

        # save
        saveCohortData(res_all,sub_ids,'%s/motion-corrected'%datadir)

