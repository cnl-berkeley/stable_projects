import os,sys
import os,sys
from numpy import *
import numpy as np
import pandas as pd

def getReqLabels(hemis,prob_tag):
    lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
    lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    mpar_labels = []
    
    labels = lpfc_labels + lpar_labels + mpar_labels
    req_labels = []
    for hemi in hemis:
        for label in labels:
            req_labels.append('%s.%s'%(hemi,label))

    return req_labels


def pruneFC(datadir,sub,req_labels):

    fc_file = '%s/sub-%s/ses-%s/seed-fc_labels_hybrid.csv'%(datadir,sub,sub[-1])
    data = pd.read_csv(fc_file)
    column_names = data.columns
    label_inds = []
    for i in arange(len(column_names)):
        if column_names[i] in req_labels:
            label_inds.append(i)
    data = array(data)
    fc_arr = data[ix_(label_inds,label_inds)]

    return fc_arr


if __name__ == '__main__':
    # combines corrected within-hemisphere lh and rh matrices and uncorrected across-hemispheric lh-rh matrices
    
    sub = sys.argv[1]
    pipeline = sys.argv[2]
    prob_tag = sys.argv[3]

    if prob_tag != 'prob':
        print('not prob')
        datadir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2seed%s/relmatch'%pipeline
        datadir_corr = '/home/weiner/shakki/NotBackedUp/spatial_autocorrelation%s/%s/labelnw'%(pipeline,sub)
    else:
        print('prob')
        datadir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2seed%s-prob_mpm/relmatch'%pipeline
        datadir_corr = '/home/weiner/shakki/NotBackedUp/spatial_autocorrelation%s/%s/labelnw_prob_mpm'%(pipeline,sub)

    req_labels = getReqLabels(['lh','rh'],prob_tag)
    data = pruneFC(datadir,sub,req_labels)

    savetxt('%s/fc-labels-both.txt'%datadir_corr,data)

    res_lh = pd.read_csv('%s/fc-labels-%s-sareg.txt'%(datadir_corr,'lh'),sep='\s+',header=None)
    res_rh = pd.read_csv('%s/fc-labels-%s-sareg.txt'%(datadir_corr,'rh'),sep='\s+',header=None)

    n_labels = int(len(req_labels)/2.0)
    data[0:n_labels,0:n_labels] = res_lh
    data[n_labels::,n_labels::] = res_rh

    savetxt('%s/fc-labels-both-sareg.txt'%datadir_corr,data)

