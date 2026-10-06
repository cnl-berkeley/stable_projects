import os,sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import pearsonr


npz = np.load('/home/weiner/shakki/scripts-new/visualization/matrix/plotting_order_man.npz')
order = npz['order']
module_colors = npz['node_colors']
labels = npz['labels']
c = npz['c']
order = order[::-1]

n_names = len(labels)
labels_ordered = np.array(labels)[order]

lpfc_labels = ['ifs','pimfs','pmfs-a','pmfs-i','pmfs-p','sfs-a','sfs-p'] + ['prts','lfms','aalf']
lpar_labels = ['slocs-v','sB','pips','mTOS','lTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
mpar_labels = []
plot_names = lpfc_labels + lpar_labels + mpar_labels

pts = ['pimfs','pmfs-a','pmfs-i','pmfs-p','aipsJ','slocs-v']

labels2 = []
for hemi in ['lh','rh']:
    for name in plot_names:
        if name in pts:
            name = '*%s'%name
        labels2.append('%s.%s'%(hemi,name))


def getOrigData(datadir,sub,label1,labels_sel2):
    fc_file = '%s/fc_seed2seed-bp/relmatch/sub-%s/ses-%s/seed-fc_labels_hybrid.csv'%(datadir,sub,sub[-1])
    data = pd.read_csv(fc_file)
    column_names = np.array(data.columns)
    idx1 = np.argwhere(column_names == label1)
    data = np.array(data)[idx1,:].squeeze()

    idx = []
    for label_sel in labels_sel2:
        if label_sel != label1:
            idx.append(np.argwhere(column_names == label_sel)[0])
    idx = np.array(idx).flatten()

    return data[idx]


def getProbData(datadir,sub,label1,labes_sel2):
    fc_file = '%s/fc_seed2seed-bp-prob_mpm/relmatch/sub-%s/ses-%s/seed-fc_labels_hybrid.csv'%(datadir,sub,sub[-1])
    data = pd.read_csv(fc_file)
    column_names = np.array(data.columns)
    idx1 = np.argwhere(column_names == label1)
    data = np.array(data)[idx1,:].squeeze()

    idx = []
    for label_sel in labels_sel2:
        if label_sel != label1:
            print(label_sel)
            idx.append(np.argwhere(column_names == label_sel)[0])
    idx = np.array(idx).flatten()

    return data[idx]


def getCorrSub(datadir,sub,labels_sel2):
    vals_man = []
    vals_prob = []
    for label1 in labels_sel2:
        try:
            vals_man.extend(getOrigData(datadir,sub,label1,labels_sel2))
            vals_prob.extend(getProbData(datadir,sub,label1,labels_sel2))
        except:
            print('error with %s'%subs[s])

    r,p = pearsonr(vals_man,vals_prob)

    return r


def doOneSulcus(datadir,subs,label1,labels_sel2):

    # collect data
    corrs = []
    for s in np.arange(n_subs):
        try:
            vals_man = getOrigData(datadir,subs[s],label1,labels_sel2)
            vals_prob = getProbData(datadir,subs[s],label1,labels_sel2)
            r,p = pearsonr(vals_man,vals_prob)
            corrs.append(r)
        except:
            print('error with %s'%subs[s])

    return corrs


if __name__ == '__main__':

    outdir = '%s/density_plots'%os.getcwd()
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p','prts','lfms','aalf']
    lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    labels = np.array(lpfc_labels + lpar_labels)
    labels_sel = []
    for hemi in ['lh','rh']:
        for label in labels:
            labels_sel.append('%s.%s'%(hemi,label))
    labels_sel = np.array(labels_sel)

    excl_sulci = ['lh.slos1','rh.pmfs_a','lh.aipsJ','rh.aipsJ']
    labels_sel2 = [value for value in labels_sel if value not in excl_sulci]
    labels_sel2 = np.array(labels_sel2)

    datadir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data'

    namelistfile = '/home/weiner/shakki/scripts-new/preprocessing/regress_motion/-bp/orig/fc-sareg-fdreg/namelist.csv'
    subs = np.loadtxt(namelistfile,dtype='str')
    n_subs = len(subs)
    print('n_subs',n_subs,subs)

    # correlation across all connections per subject
    sub_corrs = []
    for sub in subs:
        sub_corrs.append(getCorrSub(datadir,sub,labels_sel2))
    print(np.nanmean(sub_corrs),np.nanmin(sub_corrs),np.nanmax(sub_corrs))


    # calculate correlation between man and prob fingerprint per sulcus
    n_labels = len(labels_sel)
    arr_corr = np.zeros((n_labels,n_subs))*np.nan
    for l1 in np.arange(n_labels):
        label1 = labels_sel[l1]
        if label1 in labels_sel2:
            corr = doOneSulcus(datadir,subs,label1,labels_sel2)
            arr_corr[l1,:] = corr

    my_pal = module_colors[order][::-1]

    my_order = list(np.array(labels_sel)[order][::-1])

    # plot
    plt.figure()
    df = pd.DataFrame(arr_corr.flatten(),columns=['Correlation'])
    df['sulcus'] = np.array(labels_sel).repeat(n_subs)
    df['sub'] = np.tile(np.arange(n_subs),n_labels)

    ax = sns.boxplot(data=df,x="sulcus",y="Correlation",palette=my_pal,order=my_order,fliersize=2)
    ax.set_xticks(np.arange(len(labels2)),labels=my_order)
    ax.tick_params(axis='x',rotation=90)

    ax.spines['left'].set_position(('data',-1.0))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xlabel('') 

    plt.savefig('%s/tmp.png'%os.getcwd(),bbox_inches='tight')
