import os,sys
from numpy import *
from numpy.ma import masked_where
import numpy.ma
from matplotlib import *
import matplotlib.pyplot as plt
import pylab as P
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import glob
import pandas as pd
import networkx as nx
import bct

from collections import Counter
import seaborn as sns

from matplotlib.font_manager import FontProperties

npz = load('out/plotting_order_man_pTS.npz')
order = npz['order']
labels = npz['labels']
c = npz['c']
n_names = len(labels)

ciu = zeros((42,3))
for i in arange(n_names):
   ciu[i,:] = c[argwhere(order == int(i)),0,:]


def getLabels(hemis):
    label_type_groups = []
    lpfc_labels = ['ifs','pimfs*','pmfs-a*','pmfs-i*','pmfs-p*','sfs-a','sfs-p'] + ['prts*','lfms*','aalf*']
    lpar_labels = ['slocs-v*','sB','pips','mTOS','lTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ*']
    plot_names = lpfc_labels + lpar_labels

    labels = []
    label_dict = {}
    counter = 0
    for hemi in hemis:
        for name in plot_names:
            labels.append('%s.%s'%(hemi,name))
            label_dict[counter] = "%s.%s"%(hemi,name)
            counter = counter+1

    return label_dict,labels,label_type_groups


def getPos(hemis):
    # symmetric hemispheres
    a_lh = pd.read_csv('/home/weiner/shakki/scripts-new/visualization/brainnet/atlas_inflated_21labels_lh.csv',sep='\s+',header=None)
    a_rh = pd.read_csv('/home/weiner/shakki/scripts-new/visualization/brainnet/atlas_inflated_21labels_rh.csv',sep='\s+',header=None)
    ad_lh = array(a_lh)[:,0:3] 
    ad_rh = array(a_rh)[:,0:3] 
    ad = mean((ad_lh,ad_rh),axis=0)

    pos = {}
    if hemis == ['lh']:
        for i in arange(ad.shape[0]):
            pos[i] = (1-a[i][2]*150,a[i][1]*100)
    elif hemis == ['rh']:
        for i in arange(ad.shape[0]):
            pos[i] = (14000+a[i][2]*150,a[i][1]*100)
    else:         
        for i in arange(ad.shape[0]):
            pos[i] = (1-ad[i][2]*150,ad[i][1]*100)
        for i in arange(ad.shape[0]):
            pos[i+ad.shape[0]] = (14000+ad[i][2]*150,ad[i][1]*100)

    return pos


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
    if n_uni_mod > 4:
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


def do(A,thr,t_ind,label_dict,labels,label_type_groups,pos,cius,module_colors,node_sizes,savedir,savename):
    fig = P.figure(figsize=(5,4))
    a = fig.add_axes([0, 0, 1, 1])

    A_orig = A
    A[A < 0] = 0
    data_thr = bct.threshold_proportional(A,0.12,copy=True)

    for i in arange(data_thr.shape[0]):
       if sum(data_thr[i,:]) == 0:
           data_thr[i,argmax(A_orig[i,:])] = 0.001
    print(data_thr)
    A2 = data_thr

    label_dict = {}
    lpfc_labels = ['ifs','pimfs*','pmfs-a*','pmfs-i*','pmfs-p*','sfs-a','sfs-p'] + ['prts*','lfms*','aalf*']
    lpar_labels = ['slocs-v*','sB','pips','mTOS','lTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ*']
    plot_names = lpfc_labels + lpar_labels
    counter2 = 0
    for hemi in ['lh','rh']:
        for name in plot_names:
            label_dict[counter2] = "%s.%s"%(hemi,name)
            counter2 = counter2+1
    print(len(label_dict))

    # add edges
    G = nx.from_numpy_array(A2)
    G.edges(data=True)
    pos_spring = nx.kamada_kawai_layout(G,weight=None,scale=1)
    edges = nx.draw_networkx_edges(G,pos_spring,alpha=1.0,width=1,edge_color='lightgrey')
    edges.set(linestyles='dotted')

    data_thr = bct.threshold_proportional(A,0.12,copy=True)
    G = nx.from_numpy_array(data_thr)
    edges = nx.draw_networkx_edges(G,pos_spring,alpha=1.0,width=1,edge_color='lightgrey')
    edges.set(linestyles='solid')

    # add nodes
    node_colors = ciu
    node_sizes2 = array([(2**5)*abs(v) for v in node_sizes])
    node_sizes2 = node_sizes2/amax(node_sizes2)*350

    # add node borders
    nx.draw_networkx_nodes(G,pos_spring,node_color=node_colors,node_size=node_sizes2,alpha=1.0)
    nx.draw_networkx_labels(G, pos_spring, label_dict, font_size=6, font_color="black")

    ax = plt.gca()
    ax.margins(0.03)
    plt.axis("off")

    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks([])
    a.set_yticks([])

    fig.savefig('%s/%s.svg'%(savedir,savename), format="svg",bbox_inches='tight',pad_inches=0.03, dpi=600)
    plt.close()


def plotOne(tag,norm,pos,label_dict,labels,thrs):
    pipeline = '-bp'
    pipel = 'orig'

    savedir='%s/out/spring'%os.getcwd()
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    # data for edges
    datafile = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-both/regs/data/mean.csv'%(pipeline,pipel,tag,norm)
    A = loadtxt(datafile,delimiter=',')

    # data for node size
    datafile = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-both/regs/graph-nopers-bin/mean-str-auc.txt'%(pipeline,pipel,tag,norm)
    node_sizes = loadtxt(datafile,delimiter=' ')

    # get infomap modules
    datafile = '/home/weiner/shakki/scripts-new/analysis/output/infomap/%s/%s/%s-%s/mean/infomap-mean_%s_%s-fixed.txt'%(pipeline,pipel,tag,norm,norm,'both')
    cius,module_colors = getModules(datafile)
    # assign nodes without any connections 0
    res = getIsolated(A,thrs)
    cius[res == 0] = 0

    t_ind = 11
    thr = thrs[t_ind]
    savename = 'spring_%.2f'%thr

    do(A,thr,t_ind,label_dict,labels,label_type_groups,pos,cius,module_colors,node_sizes,savedir,'map-%s'%savename)


if __name__ == '__main__':

    thrs = arange(0.01,0.5,0.01)

    # node positions
    pos = getPos(['lh','rh'])

    # node names etc
    label_dict,labels,label_type_groups = getLabels(['lh','rh'])

    for tag in ['-sareg-fdreg']:
        for norm in ['nonorm']:
            plotOne(tag,norm,pos,label_dict,labels,thrs)


