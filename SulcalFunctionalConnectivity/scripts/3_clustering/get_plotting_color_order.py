import os,sys
from numpy import *
import bct
import pandas as pd
from collections import Counter
import seaborn as sns

lpfc_labels = ['ifs','pimfs','pmfs-a','pmfs-i','pmfs-p','sfs-a','sfs-p'] + ['prts','lfms','aalf']
lpar_labels = ['slocs-v','sB','pips','mTOS','lTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
mpar_labels = []
plot_names = lpfc_labels + lpar_labels + mpar_labels

labels = []
for hemi in ['lh','rh']:
    for name in plot_names:
        labels.append('%s.%s'%(hemi,name))
        
n_names = len(labels)


def getIsolated(data,thrs):
    res = zeros((shape(data)[0],len(thrs)))
    for t in arange(len(thrs)):
        data_thr = bct.threshold_proportional(data,thrs[t],copy=True)
        res[bct.strengths_und(data_thr) != 0.0,t] = 1

    return res    


def getModules(datafile):
    a = pd.read_csv(datafile,sep=',',header=None)
    a = a.astype(int)
    a = array(a)

    # modules with > 1 node
    for t in arange(shape(a)[1]):
        c = Counter(array(a[:,t]).flatten())
        for mod in c:
            if c[mod] == 1:
                a[a[:,t] == mod,t] = 0
    uni_mods = unique(a)
    n_uni_mod = len(uni_mods)

    # define colors (assuming there will be single-node modules!)
    module_colors = zeros((n_uni_mod,3))
    module_colors[0,:] = 0.9 # light grey
    module_colors[1::,:] = sns.color_palette('hls',n_uni_mod-1) # non-grey colors
    module_colors[4,:] = sns.color_palette("bright")[0]
    return a,module_colors


def getNodeColors(a,module_colors,t_ind):
    uni_mods = unique(a)
    modules = array(a)[:,t_ind]
    node_colors = zeros((len(modules),3))
    for m in arange(len(modules)):
        ind = argwhere(uni_mods == modules[m])[0]
        node_colors[m,:] = module_colors[ind]

    return node_colors


pipeline = 'orig'

# sorting and colorbar
datafile = '/home/weiner/shakki/scripts-new/analysis/output/infomap/-bp/%s/-sareg-fdreg-nonorm/mean/infomap-mean_nonorm_both-fixed.txt'%pipeline
cius,module_colors = getModules(datafile)

datafile = '/home/weiner/shakki/scripts-new/analysis/output/setup/-bp/%s/setup-sareg-fdreg/nonorm-both/regs/data/mean.csv'%pipeline
A = loadtxt(datafile,delimiter=',')

thrs = arange(0.01,0.5,0.01)
res = getIsolated(A,thrs)
cius[res == 0] = 0

# order data
t_ind = 14 # corresponding to the density of 0.15
res = array(cius)
res_vals = res[:,t_ind]

node_colors = getNodeColors(cius,module_colors,t_ind)

uni_cius = unique(res_vals)
uni_cius = [6,1,5,3,0]
order = []

c = zeros((shape(res_vals)[0],5,3))
for u in arange(len(uni_cius)):
    inds = where(res_vals == uni_cius[u])[0]
    tmp = A[inds][:,inds]
    ordered,indices,cost = bct.reorder_matrix(A[inds][:,inds])
    order.extend(inds[indices])

    wth = [4,1,3,2,0]
    c[inds,:,:] = module_colors[wth[u]]

c = c[order]

#savez('plotting_order.npz',order=order,node_colors=node_colors,labels=labels,ciu=res_vals,c=c)


### manual tuning

labels2 = array(labels)

# switch to red (1)
chlabels = ['rh.lfms','lh.pmfs-a']
for label in chlabels:
    inds = argwhere(labels2 == label)[0]
    res_vals[inds] = 1

# switch to blue (6)
chlabels = ['rh.pmfs-i']
for label in chlabels:
    inds = argwhere(labels2 == label)[0]
    res_vals[inds] = 6

# switch to green (5)
chlabels = ['rh.aalf','rh.prts']
for label in chlabels:
    inds = argwhere(labels2 == label)[0]
    res_vals[inds] = 5

# switch to cyan (8)
chlabels = ['rh.slocs-v','lh.slocs-v']
for label in chlabels:
    inds = argwhere(labels2 == label)[0]
    res_vals[inds] = 8

# switch to purple (17)
chlabels = ['lh.aipsJ','lh.lfms','lh.pmfs-i','rh.pmfs-p','rh.aipsJ']
for label in chlabels:
    inds = argwhere(labels2 == label)[0]
    res_vals[inds] = 17

# switch to lightblue (13)
chlabels = ['rh.cSTS3','lh.cSTS3']
for label in chlabels:
    inds = argwhere(labels2 == label)[0]
    res_vals[inds] = 13

# order:
# blue 6
# red 1
# purple 17
# lightblue 13
# green 5
# cyan 8
# yellow 3

#uni_cius = unique(res_vals)
uni_cius = [6,1,17,5,13,8,3]
order = []

# 6 -> 4
# 1 -> 1
# 5 -> 3
# 3 -> 2
# 0 -> 0

c = zeros((shape(res_vals)[0],5,3))
module_colors2 = zeros_like(module_colors)
print(shape(module_colors))
for u in arange(len(uni_cius)):
    inds = where(res_vals == uni_cius[u])[0]
    tmp = A[inds][:,inds]
    ordered,indices,cost = bct.reorder_matrix(A[inds][:,inds])
    order.extend(inds[indices])

    wth = [4,1,7,3,6,5,2]
    c[inds,:,:] = module_colors[wth[u]]
    print(shape(module_colors[wth[u]]))

    module_colors2[u,:] = module_colors[wth[u]]

node_colors = squeeze(c[:,0,:])
print(node_colors)
c = c[order]

savez('plotting_order_man.npz',order=order,module_colors=module_colors2,node_colors=node_colors,labels=labels,ciu=res_vals,c=c)


# add pTS indicators to sulcal names for plotting

import numpy as np

labels = ['lh.ifs','lh.pimfs*','lh.pmfs-a*','lh.pmfs-i*','lh.pmfs-p*','lh.sfs-a',\
 'lh.sfs-p','lh.prts*','lh.lfms*','lh.aalf*','lh.slocs-v*','lh.sB','lh.pips',\
 'lh.mTOS','lh.lTOS','lh.IPS-PO','lh.IPS','lh.cSTS1','lh.cSTS2','lh.cSTS3',\
 'lh.aipsJ*','rh.ifs','rh.pimfs*','rh.pmfs-a*','rh.pmfs-i*','rh.pmfs-p*',\
 'rh.sfs-a','rh.sfs-p','rh.prts*','rh.lfms*','rh.aalf*','rh.slocs-v*','rh.sB',\
 'rh.pips','rh.mTOS','rh.lTOS','rh.IPS-PO','rh.IPS','rh.cSTS1','rh.cSTS2',\
 'rh.cSTS3','rh.aipsJ*']

np.savez('%s/out/plotting_order_man_pTS.npz'%os.getcwd(),order=npz['order'],module_colors=npz['module_colors'],node_colors=npz['node_colors'],labels=labels,ciu=npz['ciu'],c=npz['c'])
