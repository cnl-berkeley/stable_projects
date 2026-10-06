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
import cmasher as cmr
from matplotlib import transforms
from collections import Counter
import seaborn as sns
import cmasher as cmr
from matplotlib.font_manager import FontProperties

fsize=9.4

# background brain outline
import matplotlib.image as mpimg
im_mni_flat = mpimg.imread('brain2.jpg')
im_mni_flat = im_mni_flat[53:-53,53:-53]


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(linspace(minval, maxval, n)))
    return new_cmap

cmap_mean = truncate_colormap(cmr.fusion_r,0.15,0.85)

def getLabels(hemis):
    lpfc_labels = ['ifs','pimfs','pmfs-a','pmfs-i','pmfs-p','sfs-a','sfs-p'] + ['prts','lfms','aalf']
    lpar_labels = ['slocs-v','sB','pips','mTOS','lTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    mpar_labels = []
    plot_names = lpfc_labels + lpar_labels + mpar_labels

    labels = []
    label_dict = {}
    counter = 0
    for hemi in hemis:
        for name in plot_names:
            labels.append('%s.%s'%(hemi,name)) 
            label_dict[counter] = name
            counter = counter+1

    return label_dict,labels


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
            pos[i] = (16000+a[i][2]*150,a[i][1]*100)
    else:         
        for i in arange(ad.shape[0]):
            pos[i] = (1-ad[i][1],140+ad[i][2])
        for i in arange(ad.shape[0]):
            pos[i+ad.shape[0]] = (1-ad[i][1],ad[i][2])

    for i in [0,21]:
        pos[2+i] = (pos[2+i][0]-16,pos[2+i][1]) # pmfs_a
        pos[7+i] = (pos[7+i][0]+5,pos[7+i][1]) # prts

        pos[14+i] = (pos[14+i][0]+10,pos[14+i][1]-7) #iTOS
        pos[13+i] = (pos[13+i][0]+5,pos[13+i][1]+3) #mTOS
        pos[11+i] = (pos[11+i][0]+5,pos[11+i][1]) #sB

        pos[18+i] = (pos[18+i][0],pos[18+i][1]-3) # cSTS2
        pos[19+i] = (pos[19+i][0]-2,pos[19+i][1]-5) # cSTS3

        pos[16+i] = (pos[16+i][0],pos[16+i][1]+7) # IPS
        pos[20+i] = (pos[20+i][0],pos[20+i][1]+5) # aipsJ

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


def getNodeColors(a):
    a = array(a)
    node_colors = zeros((42,3))
    node_colors[a == 0,:] = [0.8,0.8,0.8]
    node_colors[a == 1,:] = [178.0/255,102.0/255,255.0/255]
    node_colors[a == 2,:] = [0,204.0/255,0]
    node_colors[a == 3,:] = [255.0/255,153.0/255,51.0/255] # bw orange
    node_colors[a == 4,:] = [0.8,0.8,0.8]
    node_colors[a == 5,:] = [255.0/255,153.0/255,51.0/255] # bw orange

    return node_colors


def do(A_out,A_out2,label_dict,labels,label_type_groups,pos,node_colors,node_sizes,savedir,savename,inc_contralateral):
    fig = P.figure(figsize=(5,7))

    # upper
    a = fig.add_axes([-0.05, 0.18, 1.1, 1.1])
    tr = transforms.Affine2D().rotate_deg(2)
    a.imshow(im_mni_flat, transform=tr + a.transData,cmap=plt.cm.Greys,interpolation='bilinear',aspect='equal',alpha=0.2)
    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks([])
    a.set_yticks([])

    # lower
    a = fig.add_axes([-0.05, -0.37, 1.1, 1.1])
    a.imshow(im_mni_flat, transform=tr + a.transData,cmap=plt.cm.Greys,interpolation='bilinear',aspect='equal',alpha=0.2)
    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks([])
    a.set_yticks([])

    a = fig.add_axes([0, 0, 1, 1])

    # isolates are mispositioned so make them not isolated
    pipeline = '-bp'
    pipel = 'orig'
    norm = 'nonorm'
    tag = '-sareg-fdreg'
    datafile = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-both/regs/data/mean.csv'%(pipeline,pipel,tag,norm)
    A2 = loadtxt(datafile,delimiter=',')
    A2 = bct.threshold_proportional(A2,0.15,copy=True)
    G = nx.from_numpy_array(A2)
    G.edges(data=True)
    edges2 = nx.draw_networkx_edges(G,pos,alpha=0,width=0,edge_color='white') #edge_color,edge_cmap=cmr.fusion_r)

    # edge values
    G = nx.from_numpy_array(A_out)
    G.edges(data=True)
    edge_width=[1.2*G[u][v]['weight'] for u, v in G.edges()]
    edge_weight=[0.8*G[u][v]['weight'] for u, v in G.edges()]
    edges = nx.draw_networkx_edges(G,pos,alpha=0.9,width=2,edge_cmap=cmap_mean,edge_color=edge_weight,edge_vmin=-4.0, edge_vmax=4.0, style='solid') #edge_color,edge_cmap=cmr.fusion_r)

    # edge values
    G = nx.from_numpy_array(A_out2)
    G.edges(data=True)
    edge_width=[1.2*G[u][v]['weight'] for u, v in G.edges()]
    edge_weight=[0.8*G[u][v]['weight'] for u, v in G.edges()]
    edges2 = nx.draw_networkx_edges(G,pos,alpha=0.9,width=1.5,edge_cmap=cmap_mean,edge_color=edge_weight,edge_vmin=-4.0, edge_vmax=4.0,style=(0,(2,2))) #edge_color,edge_cmap=cmr.fusion_r)

    # nodes
    node_sizes2 = array([(2**5)*abs(v) for v in node_sizes])
    node_sizes2 = node_sizes2/amax(node_sizes2)*350

    # draw
    node_colors = getNodeColors(node_colors)
    nx.draw_networkx_nodes(G,pos,node_color='white',linewidths=0,edgecolors=None,node_size=node_sizes2,cmap=cmr.ember_r,alpha=1.0)
    nx.draw_networkx_nodes(G,pos,node_color=node_colors,linewidths=0,edgecolors=None,node_size=node_sizes2,cmap=cmr.ember_r,alpha=0.8)

    ax = plt.gca()
    ax.margins(0.03)
    plt.axis("off")

    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks([])
    a.set_yticks([])

    a = fig.add_axes([1.05,0.45,0.015,0.2])
    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks([])
    a.set_yticks([])

    a = fig.add_axes([1.05,0.75,0.015,0.2])
    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks([])
    a.set_yticks([])

    if inc_contralateral:
       savename = '%s-wcontra_sa+nosa-nolabels'%savename
    fig.savefig('%s/mod-%s.png'%(savedir,savename),bbox_inches='tight',pad_inches=0.03, dpi=200)
    plt.close()


def plotOne(tag,norm,pos,label_dict,labels,thrs,bin_tag,inc_contralateral):
    pipeline = '-bp'
    pipel = 'orig'

    savedir='/home/weiner/shakki/scripts-new/analysis/output/visualization/networkx/%s/%s/%s-%s/both-lat/%s'%(pipeline,pipel,tag,norm,bin_tag)
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    datafile = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-both/regs/graph-pers-bin/mean-str-auc.txt'%(pipeline,pipel,tag,norm)
    node_sizes = loadtxt(datafile,delimiter=' ')

    # lh.pmfs_i
    datafile = '/home/weiner/shakki/scripts-new/analysis/output/visualization/nbs/sel_regs-OLS-%s/%s/%s/depths-sorted-one-tailed.npz'%(tag,pipeline,pipel)
    npz = load(datafile)
    Aconn = npz['res_fdr_row']
    p_arr = squeeze(Aconn[:,:,1])
    t_arr = squeeze(Aconn[:,:,0])
    t_arr[p_arr > 0.05] = 0
    A_out = zeros((42,42))
    A_out[3,:] = t_arr[4,:]
    node_colors = [4]*42
    node_colors[3] = 5

    Aconn = npz['res_uncp']
    p_arr = squeeze(Aconn[:,:,1])
    t_arr = squeeze(Aconn[:,:,0])
    t_arr[p_arr > 0.05] = 0
    A_out2 = zeros((42,42))
    A_out2[2+21,:] = t_arr[3,:]
    node_colors[2+21] = 5

    savename = 'lh.pmfs_i-depth-sa-uncp'
    do(A_out,A_out2,label_dict,labels,label_type_groups,pos,node_colors,node_sizes,savedir,'map-%s'%savename,inc_contralateral)


if __name__ == '__main__':

    bin_tag = 'bin'

    thrs = arange(0.01,0.5,0.01)

    # node positions
    pos = getPos(['lh','rh'])

    # node names etc
    label_dict,labels = getLabels(['lh','rh'])

    norm = 'nonorm'

    for tag in ['-sareg-fdreg']:
        for inc_contralateral in [True,False]:
            plotOne(tag,norm,pos,label_dict,labels,thrs,bin_tag,inc_contralateral)

