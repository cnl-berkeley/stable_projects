import os,sys
import re
import numpy as np
import pandas as pd
from scipy import stats


def getData(df,tag,pipe):
    # the connectivity matrix is required to have non-nan value for all connections
    use_sub_ind = []
    data_all = []
    sub_ids = []
    sub_fds = []
    for subindex, row in df.iterrows():
        sub = row.both_id
        ses = sub[-1]
        dfile = '/home/weiner/shakki/NotBackedUp/spatial_autocorrelation%s/%s/labelnw_prob_mpm/fc-labels-both%s.txt'%(pipe,sub,tag)
        print(dfile)
        if os.path.exists(dfile):
            df_sub = pd.read_csv(dfile,sep='\s+',header=None)
            df_sub = df_sub.round(decimals=3)

            df_sub[np.isinf(df_sub)] = 0
            data_all.append(df_sub)
            sub_ids.append(sub)
            sub_fds.append(row['mFD'])

    return data_all,sub_fds,sub_ids


def getResiduals(data,correlate):
    n_subs = len(data)
    n_roi = np.shape(data[0])[1]
    x = correlate

    res_arr = np.zeros((n_subs,n_roi,n_roi))
    for i in np.arange(n_roi):
        for j in np.arange(n_roi):
            y = []
            for s in np.arange(n_subs):
                y.append(np.array(data[s])[i,j]) 
            print(y)
            ret = stats.linregress(x,y)
            for s in np.arange(n_subs):
               res_arr[s,i,j] = y[s] - ret[0]*x[s] # leave intercept

    return res_arr


def saveCohortData(data_all,sub_ids,savedir):
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    # save used subjects
    np.savetxt('%s/namelist.csv'%savedir,sub_ids,delimiter=',',header='',fmt='%s')

    # save residualized data
    all_vals = np.zeros((len(sub_ids),len(data_all[0]),len(data_all[0])))
    for i in np.arange(len(sub_ids)):
        data = np.array(data_all[i])
        savedir2 = '%s/raw/%s'%(savedir,sub_ids[i])
        if not os.path.exists(savedir2):
            os.makedirs(savedir2)
        np.savetxt('%s/fc-labels-both.csv'%savedir2,data,delimiter=',',header='',fmt='%.04f')
        all_vals[i,:,:] = np.array(data)

    # save mean residual data
    mean_data = np.squeeze(np.nanmean(all_vals,0))
    np.savetxt('%s/raw/mean_data-both.csv'%savedir,mean_data,delimiter=',',header='',fmt='%.04f')

    # save group variance residual data
    var_data = np.squeeze(np.nanvar(all_vals,0))
    np.savetxt('%s/raw/var_data-both.csv'%savedir,var_data,delimiter=',',header='',fmt='%.04f')

    # save data deviations from residual mean
    if not os.path.exists('%s/dev_from_mean'%savedir):
        os.makedirs('%s/dev_from_mean'%savedir)
    for i in np.arange(len(sub_ids)):
        data = np.array(all_vals[i])
        data = np.subtract(data,mean_data)
        savedir2 = '%s/dev_from_mean/%s'%(savedir,sub_ids[i])
        if not os.path.exists(savedir2):
            os.makedirs(savedir2)
        np.savetxt('%s/fc-labels-both.csv'%savedir2,data,delimiter=',',header='',fmt='%.04f')


if __name__ == '__main__':
    # list of potential subjects
    cohort_file = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/fmri_cohort/relmatch_parlabels3/cohort_full.csv'
    df = pd.read_csv(cohort_file)

    cohort_file2='/home/weiner/shakki/scripts-new/preprocessing/regress_motion/-bp/orig/fc-sareg/namelist.csv' # subjects in the main analysis
    subs = np.loadtxt(cohort_file2,dtype='str').tolist()
    df_filt = df[df['both_id'].isin(subs)]

    pipe = '-bp'

    pipeline = 'orig'

    savedir = '%s/prob_mpm/%s/%s'%(os.getcwd(),pipe,pipeline)
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    for tag in ['','-sareg']:
        # filter out subjects without connectivity values on required labels
        data,sub_fds,sub_ids = getData(df_filt,tag,pipe)

        # regress out effects of mfd
        res_all = getResiduals(data,sub_fds)

        # save
        saveCohortData(res_all,sub_ids,'%s/fc%s-fdreg'%(savedir,tag))

