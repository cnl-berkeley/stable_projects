import os,sys
from numpy import *
from matplotlib import *
import matplotlib.pyplot as plt
import pylab as P
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import statsmodels.stats as stats
from copy import deepcopy

import cmasher as cmr

npz = load('out/plotting_order_man_pTS.npz')
order = npz['order']
module_colors = npz['node_colors']
labels = npz['labels']
c = npz['c']

n_names = len(labels)
labels_ordered = array(labels)[order]

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(linspace(minval, maxval, n)))
    return new_cmap

cmap_mean = truncate_colormap(cmr.fusion_r,0.15,0.85)

ylabels = ['Degree','Betweenness','Participation']

def saveFig(vals_all,title,savedir,savename):
    fig = P.figure(figsize=(12,12))
    
    a = fig.add_axes([0.25, 0.25, 0.46, 0.05])
    data_t = vals_all[:,:,0]
    data_t = data_t[:,order]
    vmax = ceil(nanmax(abs(data_t)))
    im = a.imshow(data_t, extent=[0,1,0,1], aspect='auto', interpolation='nearest', cmap=cmap_mean, origin='upper',vmin=-1.0*vmax,vmax=vmax,alpha=1.0)
    lims = im.get_clim()

    x = a.get_xlim()
    y = a.get_ylim()

    data_p = vals_all[:,:,1]
    data_p[isnan(data_p)] = 100
    data_p = data_p[:,order]
    print(shape(data_p))

    n_names2 = shape(vals_all)[0]
    xd = 1.0*(x[1]-x[0])/n_names
    yd = 1.0*(y[1]-y[0])/n_names2

    for i in range(0,n_names2,1):
        for j in range(0,n_names,1):
            if data_p[i][j] < 0.05:
                if data_t[i][j] > 0:
                    a.plot(x[0]+0.5*xd+xd*j,y[1]-0.5*yd-yd*i,marker='o',markersize=4,color='black')
                else:
                    a.plot(x[0]+0.5*xd+xd*j,y[1]-0.5*yd-yd*i,marker='o',markersize=4,color='black')

    # grid lines
    y_ticks = []
    for i in range(0,n_names2,1):
        y_ticks.append(y[1]-0.5*yd-yd*i)
    x_ticks = []
    for i in range(0,n_names,1):
        x_ticks.append(x[1]-0.5*xd-xd*i)

    a.spines['left'].set_linewidth(0.5)
    a.spines['right'].set_linewidth(0)
    a.spines['top'].set_linewidth(0)
    a.spines['bottom'].set_linewidth(0.5)

    a.set_yticks(y_ticks)
    a.xaxis.tick_bottom()
    a.set_yticklabels(ylabels)
    a.set_xticks(x_ticks)
    a.set_xticklabels(labels_ordered[::-1])
    a.xaxis.set_label_position('bottom')
    plt.xticks(rotation = 90)
    a.tick_params(axis='x', length=13)
    a.tick_params(axis='y', length=3)

    a = fig.add_axes([0.25+0.0005, 0.25-0.0125, 0.46, 0.01]) # x orientation
    a.imshow(swapaxes(c,0,1), aspect='auto', interpolation='nearest')
    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks([])
    a.set_yticks([])

    a = fig.add_axes([0.72,0.25,0.008,0.05])
    cbar = plt.colorbar(im,cax=a)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=0)
    cbar.set_ticks([lims[0],0,lims[1]])
    cbar.ax.set_yticklabels(['%.1f'%lims[0],'0','%.1f'%lims[1]])
    cbar.set_label('t', labelpad=-10)

    if not os.path.exists(savedir):
        os.makedirs(savedir)

    fig.savefig('%s/%s.svg'%(savedir,savename), format="svg",bbox_inches='tight',pad_inches=0.03, dpi=600)
#    fig.savefig('%s/%s.png'%(savedir,savename),bbox_inches='tight',pad_inches=0.03, dpi=800)


def saveFigs(vals_all,title,savedir,savename):
    saveFig(vals_all,'',savedir,'%s_fdrrow'%savename)


if __name__ == '__main__':

    tag = '-sareg-fdreg'
    #tag = '-fdreg'

    savedir = '%s/out/centrality'%os.getcwd()
    n_labels = 42

    vals_all = zeros((3,n_labels,2))
    metrics = ['str','bw','pc']
    for m in arange(len(metrics)):
        vals_all[m,:,0] = array(pd.read_csv('/home/weiner/shakki/scripts-new/analysis/output/setup/-bp/orig/setup%s/nonorm-both/regs/graph-nopers-bin-tt/tvsmean-node_t-%s-auc.txt'%(tag,metrics[m]),header=None)).flatten()
        vals_all[m,:,1] = array(pd.read_csv('/home/weiner/shakki/scripts-new/analysis/output/setup/-bp/orig/setup%s/nonorm-both/regs/graph-nopers-bin-tt/tvsmean-node_%s-auc_fdr.txt'%(tag,metrics[m]),header=None)).flatten()

    savename = 'centrality'
    saveFigs(vals_all,'',savedir,savename)

