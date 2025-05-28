import os,sys
from numpy import *
from matplotlib import *
import matplotlib.pyplot as plt
import pylab as P
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
import bct
import pandas as pd
import cmasher as cmr
import colorcet as cc

import warnings
warnings.filterwarnings("ignore")

npz = load('%s/out/plotting_order_man_pTS.npz'%os.getcwd())
order = npz['order']
module_colors = npz['node_colors']
labels = npz['labels']
c = npz['c']
n_names = len(labels)
labels_ordered = array(labels)[order]


def plotScalar(fig,file1,lims,row_ind,col_ind,mymap,cb_pars):
    font0 = FontProperties(size=4)
    mx=0.2

    print('*file1',file1)
    try:
        data = loadtxt(file1,delimiter=',')
    except:
        data = loadtxt(file1,delimiter=' ')
    data = array(data)
    data = data[order][:,order]

    data[data <= 0] = 0.00001
    data[data == 1] = 0.99999

    a = fig.add_axes([0.25, 0.25, 0.46, 0.46])

    data[eye(shape(data)[0]) == 1] = 0

    data = tril(data)
    data[data == 0] = nan
    im = a.imshow(data, origin='upper',extent=[0,1,0,1], aspect='auto', interpolation='nearest', cmap=mymap, vmin=lims[0], vmax=lims[1])

    x = a.get_xlim()
    y = a.get_ylim()

    xd = 1.0*(x[1]-x[0])/n_names
    yd = 1.0*(y[1]-y[0])/n_names
    
    # grid lines
    linemarks = arange(0.5,n_names+0.5,1)
    y_ticks = []
    for i in range(0,n_names,1):
        y_ticks.append(y[1]-0.5*yd-yd*i)

    a.spines['left'].set_linewidth(0.5)
    a.spines['right'].set_linewidth(0)
    a.spines['top'].set_linewidth(0)
    a.spines['bottom'].set_linewidth(0.5)

    a.set_yticks(y_ticks)
    a.xaxis.tick_bottom()
    a.set_yticklabels(labels_ordered)
    a.set_xticks(y_ticks)
    a.set_xticklabels(labels_ordered[::-1])
    a.xaxis.set_label_position('bottom')
    plt.xticks(rotation = 90)
    a.tick_params(axis='both', length=13)
    
    a = fig.add_axes([0.72,0.25,0.01,0.13])
    cbar = plt.colorbar(im,cax=a)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=0)
    if cb_pars == None:
        cbar.set_ticks([lims[0],0,lims[1]])
        cbar.ax.set_yticklabels(['%.2f'%lims[0],'0','%.2f'%lims[1]]) 
    else:
        cbar.set_ticks(cb_pars[0])
        cbar.ax.set_yticklabels(cb_pars[1])
        cbar.set_label(cb_pars[2], labelpad=-10) # -15

    a = fig.add_axes([0.25-0.012,0.25-0.001,0.01,0.46]) # y orientation
    a.imshow(c, interpolation='nearest', aspect='auto')
    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks([])
    a.set_yticks([])

    a = fig.add_axes([0.25+0.0005, 0.25-0.0125, 0.46, 0.01]) # x orientation
    a.imshow(swapaxes(c,0,1), interpolation='nearest', aspect='auto')
    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks([])
    a.set_yticks([])

    return fig,im


def doOne(files,savepath,thrs,cmap,cb_pars=None):
    plt.rcParams['axes.facecolor']='white'
    plt.rcParams['savefig.facecolor']='white'
    plt.rcParams.update({'mathtext.default':'regular' })

    fig = P.figure(figsize=(12,12))
    fig,im = plotScalar(fig,files[0],thrs,0,0,cmap,cb_pars)
 
    print('* saving')
    fig.savefig('%s.svg'%savepath,format="svg",bbox_inches='tight',pad_inches=0.03, dpi=600)
    plt.close()


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(linspace(minval, maxval, n)))
    return new_cmap


if __name__ == '__main__':

    pipe = '-bp'

    pipeline = 'orig'

    savedir='%s/out/matrix'%os.getcwd()
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    cmap_var = cmr.lavender
    cmap_var2 = truncate_colormap(cmr.rainforest,0,0.95)

    cmap_mean = truncate_colormap(cmr.fusion_r,0.15,0.85)

    tags = ['-fdreg','-sareg-fdreg']
    for tag in tags:
        for norm in ['nonorm']:
            # mean
            files = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-both/regs/data/mean.csv'%(pipe,pipeline,tag,norm)
            savepath = '%s/mean%s-%s_both-nonorm15'%(savedir,tag,norm)
            doOne([files],savepath,[-0.7,0.7],cmap_mean,[[-0.7,0.7],['-0.7','0.7'],'Correlation (mean)'])

            # std
            files = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-both/regs/data/std.csv'%(pipe,pipeline,tag,norm)
            savepath = '%s/std%s-%s_both-nonorm15'%(savedir,tag,norm)
            #doOne([files],savepath,[0.15,0.3],cmap_var,[[0.15,0.3],['0.15','0.3'],'Correlation (SD)'])

            # co-clustering
            files = '/home/weiner/shakki/scripts-new/analysis/output/infomap/%s/%s/%s-%s/mean/infomap-%s_both-coiperc4-man.txt'%(pipe,pipeline,tag,norm,norm)
            savepath = '%s/coclustering-%s-%s'%(savedir,tag,norm)
            doOne([files],savepath,[0,0.8],cmap_var2,[[0,0.8],['0','80'],'% participants'])

    # mean distance
    files = '/home/weiner/shakki/scripts-new/analysis/output/setup/mean_adj.csv'
    data = loadtxt(files,delimiter=',')
    data[data == 0] = 100.0
    data[eye(shape(data)[0]) == 1] = 0
    savepath = '%s/adj'%savedir
    #doOne([files],savepath,[0,100],[[0,10,100],['0','1','20'],'Distance (cm)'])
    
