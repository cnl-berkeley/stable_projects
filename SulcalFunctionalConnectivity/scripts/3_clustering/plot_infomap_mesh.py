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

from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from collections import Counter
import seaborn as sns

from matplotlib.font_manager import FontProperties

def getLabels(hemis):
    label_anat_groups = []
    label_type_groups = []
    lpfc_labels = ['ifs','pimfs','pmfs-a','pmfs-i','pmfs-p','sfs-a','sfs-p'] + ['prts','lfms','aalf']
    lpar_labels = ['slocs-v','sB','pips','mTOS','lTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    mpar_labels = []
    plot_names = lpfc_labels + lpar_labels + mpar_labels

    labels = []
    label_dict = {}
    counter = 0
    for hemi in hemis:
        for name in plot_names:
            label_dict[counter] = "%s.%s"%(hemi,name)
            labels.append("%s.%s"%(hemi,name))
            counter = counter+1

    return label_dict,labels,label_type_groups



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


def plotInfomaps(res,t_ind,labels,module_colors,savedir,savename):
    res = array(res)

    # order data
    order = argsort(res[:,t_ind])
    res_ordered = res[order,:]
    y_ordered = array(labels)[order]

    # plot
    plt.rcParams['axes.facecolor']='white'
    plt.rcParams['savefig.facecolor']='white'
    fig = P.figure(figsize=(12,12))
    a = fig.add_axes([0.25, 0.25, 0.46, 0.46])

    uni_mods = unique(res)
    c = zeros((shape(res)[0],shape(res)[1],3))
    for m in arange(len(uni_mods)):
        c[res_ordered == uni_mods[m],:] = module_colors[m]
    im = a.imshow(c,interpolation='nearest',aspect=2)

    a.spines['left'].set_linewidth(0.5)
    a.spines['right'].set_linewidth(0)
    a.spines['top'].set_linewidth(0)
    a.spines['bottom'].set_linewidth(0.5)
    thrs = arange(0.01,0.5,0.01)

    def myformat(x, pos):
        y = 0.01*x
        return y

    a.xaxis.set_major_locator(MultipleLocator(10))
    a.xaxis.set_major_formatter(myformat)
    a.xaxis.set_minor_locator(MultipleLocator(1))

    a.set_yticks(arange(len(labels)))
    a.set_yticklabels(y_ordered)
    a.tick_params(axis='both', which='major', length=5)
    a.set_xlabel('Density')

    fig.savefig('%s/infomap-%s.svg'%(savedir,savename), format="svg",bbox_inches='tight',pad_inches=0.03, dpi=600)
    plt.close()


if __name__ == '__main__':

    pipe = '-bp'
    pipeline = 'orig'
    tag = '-sareg-fdreg'
    norm = 'nonorm'

    npz = load('out/plotting_order_man_pTS.npz')
    labels = npz['labels']

    thrs = arange(0.01,0.5,0.01)
    t_ind = 11

    # data for edges
    datafile = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-both/regs/data/mean.csv'%(pipe,pipeline,tag,norm)
    A = loadtxt(datafile,delimiter=',')
    res = getIsolated(A,thrs)

    # data on infomap
    datafile = '/home/weiner/shakki/scripts-new/analysis/output/infomap/%s/%s/%s-%s/mean/infomap-mean_%s_%s-fixed.txt'%(pipe,pipeline,tag,norm,norm,'both')
    cius,module_colors = getModules(datafile)
    cius[res == 0] = 0

    outdir = '%s/out/infomap'%os.getcwd()
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    plotInfomaps(cius,t_ind,labels,module_colors,outdir,'infomap3-%s_%0.2f'%(norm,thrs[t_ind]))
