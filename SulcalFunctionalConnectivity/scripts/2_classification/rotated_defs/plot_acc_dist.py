import os,sys
import numpy as np


def getLabelInds(req_labels):
    # wanted label names
    lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
    lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    mpar_labels = []    
    all_labels = lpfc_labels + lpar_labels + mpar_labels

    full_names = []
    for hemi in ['lh','rh']:
        for label in all_labels:
            full_names.append('%s.%s'%(hemi,label))
    full_names = np.array(full_names)

    # corresponding indices in the seed2annot output
    label_inds = []
    for label in req_labels:
        label_inds.append(int(np.argwhere(full_names == label)[0]))

    return label_inds


if __name__ == '__main__':

    conf_mats = []
    p_mats = []
    acc_mats = []
    for i in np.arange(1000):
        npzfile = '/home/weiner/shakki/scripts-github/analysis/classification/chance_level_analysis/perm_results/sulcal-classification-nosa-pair-scaled-%d.npz'%(i+1)
        npz = np.load(npzfile,allow_pickle=True)
        
        conf_mat = npz['conf_mat']
        total_acc = npz['total_acc']
        total_p = npz['total_p']
        total_acc[np.eye(42) > 0] = np.nan

        conf_mats.append(conf_mat)
        p_mats.append(total_p)
        acc_mats.append(total_acc)

    conf_mat = np.mean(np.stack(conf_mats),axis=0)
    total_p = np.mean(np.stack(p_mats),axis=0)
    total_acc = np.mean(np.stack(acc_mats),axis=0)
    total_acc[np.eye(np.shape(total_acc)[0]) == 1] = np.nan

    ###

    print('ORIGINAL')

    npzfile = '/home/weiner/shakki/scripts-new/analysis/group-analysis/#ml/sklearn/sulcal-classification-nosa-pair-scaled.npz'
    npz = np.load(npzfile,allow_pickle=True)
    conf_mat = npz['conf_mat']
    total_acc = npz['total_acc']
    total_acc[np.eye(np.shape(total_acc)[0]) == 1] = np.nan
    total_p = npz['total_p']


    # plot
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    plt.figure()
    orig_vec = total_acc.flatten()
    rot_vec = np.stack(acc_mats,axis=0).flatten()
    arr = np.concatenate((rot_vec,orig_vec))
    df = pd.DataFrame(arr,columns=['Accuracy'])
    df['measure'] = np.array(['Rotated labels']*len(rot_vec) + ['Sulcal labels']*len(orig_vec))
    g = sns.displot(df,x='Accuracy',hue='measure',stat='density',common_norm=False,bins=50,height=4,aspect=1)
    ax = plt.gca()
    ax.set_xlim(0.2,1.01)
    plt.xlabel("Classification accuracy")
    plt.ylabel("Pairwise sulcal classifications\nfor all sulci and participants") 
    g._legend.set_title('')
    plt.savefig('%s/hist_acc-polished.png'%os.getcwd(),dpi=900)
