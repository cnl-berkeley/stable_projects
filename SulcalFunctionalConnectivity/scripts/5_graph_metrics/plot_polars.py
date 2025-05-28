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

import matplotlib.pyplot as plt

import cmasher as cmr

npz = load('out/plotting_order_man_pTS.npz')
order = npz['order']
module_colors = npz['node_colors']
labels = npz['labels']
order = order[::-1]
labels2 = list(labels)

def plotOne(vals,module_colors2,module_hatches,savedir,savename,draw_labels=False):
    fig = P.figure(figsize=(12,12))

    # Compute pie slices
    N = shape(vals)[1]
    theta = linspace(0.0, 2*pi, N, endpoint=False)
    width = 2*pi/N
    colors = module_colors2[order]

    hatches = array(module_hatches)[order]

    height = nanmean(vals,0)[order]
    yerr = nanstd(vals,0)[order]

    ax = plt.subplot(projection='polar')
    container = ax.bar(theta, height, width=width, yerr=yerr, color=colors, hatch=hatches,alpha=0.9, zorder=1) #, error_kw=dict(ecolor='grey', lw=1, capsize=5, capthick=1))

    connector, caplines, (vertical_lines,) = container.errorbar.lines
    vertical_lines.set_color(colors)
    vertical_lines.set_linewidth(4)
    vertical_lines.set_alpha(1)
    vertical_lines.set_zorder(0)

    # hide grid
    ax.spines['polar'].set_visible(False)
    ax.yaxis.grid(linewidth=1,linestyle='dotted') #linestyle='.', linewidth=1)

    fsize = 34

    # mark labels
    ax.set_xticks(theta)
    ax.set_xticklabels(labels[order],fontsize=fsize,zorder=10)

    ax.set_ylim([-0.2,0.8])

    ax.set_rgrids([0, 0.2, 0.4, 0.6],labels=['','','',''],zorder=0)

    plt.gcf().canvas.draw()    
    if draw_labels:
        angles = linspace(0,2*pi,len(ax.get_xticklabels())+1)
        angles[cos(angles) < 0] = angles[cos(angles) < 0] + pi
        angles = rad2deg(angles)
        labels_new = []
        for label, angle in zip(ax.get_xticklabels(), angles):
            x,y = label.get_position()
            lab = ax.text(x,y, label.get_text(), transform=label.get_transform(), # DRAW LABELS!!
                          ha=label.get_ha(), va=label.get_va(), fontsize=fsize)
            lab.set_rotation(angle)
            labels_new.append(lab)
    ax.set_xticklabels([])

    ax.fill_between(linspace(0,2*pi,100), -0.001, 0.001, color='black', zorder=10) # <-- Added here

    fig.savefig('%s/%s.svg'%(savedir,savename), format="svg",bbox_inches='tight',pad_inches=0.03, dpi=600)


if __name__ == '__main__':

    pipe = '-bp'

    tag = '-sareg-fdreg'  #'-fdreg'

    pipeline = 'orig'

    savedir = '%s/out/polars/%s'%(os.getcwd(),tag)
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    npz = load('/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/nonorm-both/regs/data/allsub.npz'%(pipe,pipeline,tag))
    data = npz['data']

    for label in ['pimfs*','pmfs-i*','pmfs-a*']:
        for hemi in ['lh','rh']:
            print(labels2,label)
            ind = labels2.index('%s.%s'%(hemi,label))
            vals = squeeze(data[:,ind,:])
            vals[:,ind] = 0
            module_hatches=[None]*len(module_colors)
            plotOne(vals,module_colors,module_hatches,savedir,'%s.%s'%(hemi,label),True)


