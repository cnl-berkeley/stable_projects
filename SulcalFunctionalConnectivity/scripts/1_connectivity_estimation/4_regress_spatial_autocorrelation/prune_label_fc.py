import os,sys
import numpy as np
import pandas as pd


def getReqLabels(hemis,labels):
    req_labels = []
    for hemi in hemis:
        for label in labels:
            req_labels.append('%s.%s'%(hemi,label))

    return req_labels


def pruneFC(datadir,sub,hemi,req_labels,savedir):
    fc_file = '%s/sub-%s/ses-%s/seed-fc_labels_hybrid.csv'%(datadir,sub,sub[-1])
    print(fc_file)
    data = pd.read_csv(fc_file)
    column_names = data.columns
    label_inds = []
    for i in np.arange(len(column_names)):
        if column_names[i] in req_labels:
            label_inds.append(i)
    data = np.array(data)
    fc_arr = data[np.ix_(label_inds,label_inds)]
    print(fc_arr[0:5,0])
    np.savetxt('%s/fc-labels-%s.txt'%(savedir,hemi),fc_arr)


if __name__ == '__main__':

    sub = sys.argv[1]
    hemi = sys.argv[2]
    pipeline = sys.argv[3]
    prob_tag = sys.argv[4]

    if prob_tag != 'prob':
        print('not prob')
        datadir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2seed%s/relmatch'%pipeline
        outdir = '/home/weiner/shakki/NotBackedUp/spatial_autocorrelation%s/%s/labelnw2'%(pipeline,sub) 

        lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
        lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
        mpar_labels = []
        labels = lpfc_labels + lpar_labels + mpar_labels

    else:
        print('prob')
        datadir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2seed%s-prob_mpm/relmatch'%pipeline
        outdir = '/home/weiner/shakki/NotBackedUp/spatial_autocorrelation%s/%s/labelnw_prob_mpm'%(pipeline,sub) 

        lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
        lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
        mpar_labels = []
        labels = lpfc_labels + lpar_labels + mpar_labels
      
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    req_labels = getReqLabels([hemi],labels)
    pruneFC(datadir,sub,hemi,req_labels,outdir)

