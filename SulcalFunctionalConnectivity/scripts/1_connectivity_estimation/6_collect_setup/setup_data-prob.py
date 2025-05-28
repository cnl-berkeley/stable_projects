import os,sys
from numpy import *
import pandas as pd
import numpy as np
import bct
from copy import deepcopy

def getData(df,correlate_names,pipe,pipeline,tag):
    sub_ids = []
    data_all = []
    correlates_all = []
    for subindex, row in df.iterrows():
        dfile = '/home/weiner/shakki/scripts-github/analysis/network_generation/5_regress_motion/prob_mpm/%s/%s/fc%s/raw/%s/fc-labels-both.csv'%(pipe,pipeline,tag,row.both_id)
        if os.path.exists(dfile):
            df_sub = pd.read_csv(dfile,sep=',',header=None)
            df_sub = df_sub.round(decimals=3)

            correlates = []
            for correlate_name in correlate_names:
                correlates.append(row[correlate_name])

            df_sub[isinf(df_sub)] = 0
            data_all.append(df_sub)
            sub_ids.append(row.both_id)

            correlates_all.append(correlates)
        else:
            print('does not exist',dfile)


    correlates_all = array(correlates_all).T 
    return data_all,sub_ids,correlates_all


def selHemis(data_orig,hemis):
    data = deepcopy(data_orig)
    print(shape(data))
    n_subs = len(data)
    for s in arange(n_subs):
         if hemis == 'both':
             A = data[s]
         elif hemis == 'lh':
             nn = int(shape(data[s])[0]/2)
             A = array(data[s])[0:nn,0:nn]
         else:
             nn = int(shape(data[s])[0]/2)
             A = array(data[s])[nn:(2*nn),nn:(2*nn)]
         data[s] = A

    return data


def writeSetup(data,sub_ids,correlates,correlate_names,outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for c in arange(len(correlates)):
        savetxt('%s/%s.csv'%(outdir,correlate_names[c]),correlates[c],delimiter=',',header='',fmt='%.04f')

    savetxt('%s/namelist.csv'%outdir,sub_ids,delimiter=',',header='',fmt='%s')

    outdir = '%s/data'%outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for i in arange(len(sub_ids)):
        savetxt('%s/1%04d.csv'%(outdir,i),data[i],delimiter=',',header='',fmt='%.04f')

    # save mean
    mean_mat = mean(data,0) 
    savetxt('%s/mean.csv'%outdir,mean_mat,delimiter=',',header='',fmt='%.04f')

    # save std
    var_mat = std(data,0) 
    savetxt('%s/std.csv'%outdir,var_mat,delimiter=',',header='',fmt='%.04f')

    # individual data matrices
    savez('%s/allsub.npz'%(outdir),data=data)
    print(outdir)


if __name__ == '__main__':
    pipe = '-bp'

    pipeline = 'orig'

    cohort_file = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/fmri_cohort/relmatch_parlabels3check/cohort_full.csv'
    df = pd.read_csv(cohort_file)

    cohort_file2 = '/home/weiner/shakki/scripts-new/preprocessing/regress_motion/-bp/orig/fc-sareg/namelist.csv' # subjects in the main analysis
    subs = np.loadtxt(cohort_file2,dtype='str').tolist()
    df_filt = df[df['both_id'].isin(subs)]
    print(df_filt)

    df.rename(columns = {'sulcal_depth_max_lh.painfs_combined':'sulcal_depth_max_lh.pimfs_any',\
    'sulcal_depth_max_rh.painfs_combined':'sulcal_depth_max_rh.pimfs_any'}, inplace = True)

    # select samples with data and correlates
    tags = ['-sareg-fdreg','-fdreg']
    for tag in tags:
        correlate_names = ['ScanAge','mFD']
        for hemi in ['lh','rh']:
            for sulcus in ['pimfs_any','pmfs_a','pmfs_i','pmfs_p','aipsJ','slos1']:
                correlate_names.append('sulcal_depth_max_%s.%s'%(hemi,sulcus))        
        data_all,sub_ids,correlates = getData(df,correlate_names,pipe,pipeline,tag)
        print(len(sub_ids))

        savedir = '/home/weiner/shakki/scripts-github/analysis/network_generation/5_regress_motion/prob_mpm/%s/%s'%(pipe,pipeline)

        for hemis in ['both','lh','rh']:
            data_hemis = selHemis(data_all,hemis)
            writeSetup(data_hemis,sub_ids,correlates,correlate_names,'%s/setup%s/nonorm-%s/regs'%(savedir,tag,hemis))

