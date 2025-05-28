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

from statsmodels.stats.multitest import multipletests

from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()

import cmasher as cmr

npz = load('/home/weiner/shakki/scripts-new/visualization/matrix/plotting_order_man.npz')
order = npz['order']
module_colors = npz['node_colors']
labels = npz['labels']
c = npz['c']


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(linspace(minval, maxval, n)))
    return new_cmap
cmap_mean = truncate_colormap(cmr.fusion_r,0.15,0.85)

n_names = len(labels)
labels_ordered = array(labels)[order]

lpfc_labels = ['ifs','pimfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
mpar_labels = []
plot_names = lpfc_labels + lpar_labels + mpar_labels

labels2 = []
for hemi in ['lh','rh']:
    for name in plot_names:
        labels2.append('%s.%s'%(hemi,name))

ylabels = []
for name in ['pmfs_p','aipsJ','slos1']: #'pimfs','pmfs_a','pmfs_i'
    for hemi in ['lh','rh']:
        ylabels.append('%s.%s'%(hemi,name))


import statsmodels.api as sm

def regression(df,x_names,y_names,tail):
    x = df[x_names]
    y = df[y_names]

    x = sm.add_constant(x)
    try:
        model = sm.OLS(y,x)
        res = model.fit()
        
        t = res.tvalues[x_names[0]]
        p = res.pvalues[x_names[0]]
        # onetail
        if tail == 'one-tailed':
            if t > 0:
                p = 0.5*p
            else:
                p = 1.0  
    except:
        t = nan
        p = nan

    return t,p


def multiCorr(p_vals):
    res = multipletests(p_vals, alpha=0.05, method='fdr_bh') 
    return res[1]


def saveFig(vals_all,sulcus_plot_names,savedir,savename):
    fig = P.figure(figsize=(12,12))
    
    data = vals_all[:,:,0]
    
    lims = [-4,4]
    data[data < lims[0]] = lims[0]+0.0001
    data[data > lims[1]] = lims[1]-0.0001

    a = fig.add_axes([0.25, 0.25, 0.46, 0.1])
    data[data == 0] = nan
    data = data[:,order]
    im = a.imshow(data, extent=[0,1,0,1], aspect='auto', interpolation='nearest', cmap=cmap_mean, origin='upper',vmin=-4,vmax=4,alpha=0.9)
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
            	a.plot(x[0]+0.5*xd+xd*j,y[1]-0.5*yd-yd*i,marker='o',markersize=4,color='k')

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

    print('x',x_ticks,labels_ordered[::-1])
    print('y',y_ticks,ylabels)

    a.set_yticks(y_ticks)
    a.xaxis.tick_bottom()
    ylabels2 = sulcus_plot_names
    a.set_yticklabels(ylabels2,fontsize=10)
    a.set_xticks(x_ticks)
    a.set_xticklabels(labels_ordered[::-1],fontsize=10)
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

    a = fig.add_axes([0.72,0.25,0.01,0.13])
    cbar = plt.colorbar(im,cax=a)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=0)
    cbar.set_ticks([lims[0],0,lims[1]])
    cbar.ax.set_yticklabels(['%.1f'%lims[0],'0','%.1f'%lims[1]],fontsize=10)
    cbar.set_label('$t$ value', labelpad=-10)

    if not os.path.exists(savedir):
        os.makedirs(savedir)

    fig.savefig('%s/%s-fix.png'%(savedir,savename),bbox_inches='tight',pad_inches=0.03, dpi=300)


def plotScatter(df_label,t,p,col1,col2,name1,name2,ylabel,xlabel,savedir,label1,label2):

    fig = P.figure(figsize=(5,5))

    a = fig.add_axes([0.25, 0.25, 0.46, 0.46])

    if p < 0.05:
        print(label1,label2,t,p)

    print(df_label)
    print(col1,col2)
    sns.regplot(x=col1,y=col2,data=df_label,y_partial='ScanAge',ax=a,scatter_kws={'s':15})

    a.text(0.65,0.14,'p = %.3f'%p, transform=a.transAxes)
    a.text(0.65,0.07,'t = %.2f'%t, transform=a.transAxes)
    a.set_xlabel(xlabel, size = 11)
    a.set_ylabel(ylabel, size = 11)
    a.spines['left'].set_linewidth(0.5)
    a.spines['right'].set_linewidth(0)
    a.spines['top'].set_linewidth(0)
    a.spines['bottom'].set_linewidth(0.5)

    plt.savefig('%s/scatter-%s-%s.png'%(savedir,name1,name2),bbox_inches='tight')
    plt.clf()


def procSulci(sulci,sulcus_plot_names):

    pipe = '-bp'

    tag = '-sareg-fdreg' #'-fdreg'

    pipeline = 'orig'

    savedir = '/home/weiner/shakki/scripts-new/analysis/output/visualization/nbs/sel_regs-OLS-%s/%s/%s'%(tag,pipe,pipeline)
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    tail = 'one-tailed'

    npz = load('/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/nonorm-both/regs/data/allsub.npz'%(pipe,pipeline,tag))
    data = npz['data']
    n_subs = shape(data)[0]
    n_labels = shape(data)[1]

    reg = 'ScanAge'
    age_vals = array(pd.read_csv('/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/nonorm-both/regs/%s.csv'%(pipe,pipeline,tag,reg),header=None)).flatten()

    res_all = zeros((len(sulci)*2,n_labels,4))*nan
    counter = 0
    for label in sulci: # depth of which sulcus
        for hemi in ['lh','rh']:
            ind = labels2.index('%s.%s'%(hemi,label))
            vals = squeeze(data[:,ind,:])
            reg = 'sulcal_depth_max_%s.%s'%(hemi,label)
            reg_vals = array(pd.read_csv('/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/nonorm-both/regs/%s.csv'%(pipe,pipeline,tag,reg),header=None)).flatten()

            for l in arange(n_labels): # all sulci
                df = pd.DataFrame([vals[:,l],reg_vals,age_vals])
                df = df.T
                df.columns = ['fc',reg,'ScanAge']
    
                df[reg] = df[reg].astype(float)
                df['ScanAge'] = df['ScanAge'].astype(float)
                df[reg] = scaler.fit_transform(df[reg].values.reshape(-1,1))
                df['ScanAge'] = scaler.fit_transform(df['ScanAge'].values.reshape(-1,1))

                res_all[counter,l,0:2] = regression(df,[reg,'ScanAge'],['fc'],tail)

            p_vals = squeeze(res_all[counter,:,1]).flatten() # correct across row

            res_all[counter,:,2] = multiCorr(p_vals)

            counter = counter+1


    # plot scatterplots
    ind = labels2.index('lh.pmfs_i')
    vals = squeeze(data[:,ind,:])
    reg = 'sulcal_depth_max_lh.pmfs_i'
    reg_vals = array(pd.read_csv('/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/nonorm-both/regs/%s.csv'%(pipe,pipeline,tag,reg),header=None)).flatten()
    for l in arange(n_labels):
        if labels2[l] == 'lh.sfs_p':
            tlabel = 'lh.sfs-p'
        elif labels2[l] == 'lh.pmfs_p':
            tlabel = 'lh.pmfs-p'
        elif labels2[l] == 'lh.iTOS':
            tlabel = 'lh.lTOS'
        elif labels2[l] == 'rh.iTOS':
            tlabel = 'rh.lTOS'
        else:
            tlabel = labels2[l]
        print('**')
        df = pd.DataFrame([vals[:,l],reg_vals,age_vals])
        df = df.T
        df.columns = ['fc',reg,'ScanAge']
        df[reg] = df[reg].astype(float)
        df['fc'] = df['fc'].astype(float)
        df['ScanAge'] = df['ScanAge'].astype(float)

        ylabel = 'FC (%s, %s)'%('lh.pmfs-i',tlabel)
        xlabel = 'Depth (lh.pmfs-i)'
        t = res_all[4,l,0]
        p = res_all[4,l,2]
        if p < 0.05:
            plotScatter(df,t,p,'fc',reg,'scatter-marg-%s-%s'%(labels2[l],reg),reg,xlabel,ylabel,savedir,'lh.pmfs_i',labels2[l])


    # make colormesh plot
    res_all = res_all[:,:,[0,2]]
    savename = 'depths-sorted-%s-wage-%s'%(sulci[0],tail)
    saveFig(res_all,sulcus_plot_names,savedir,savename)


if __name__ == '__main__':

    sulcus_plot_names = ['lh.pimfs','rh.pimfs','lh.pmfs-a','rh.pmfs-a','lh.pmfs-i','rh.pmfs-i'] # only for plotting
    procSulci(['pimfs_any','pmfs_a','pmfs_i'],sulcus_plot_names)

    sulcus_plot_names = ['lh.pmfs-p','rh.pmfs-p','lh.aipsJ','rh.aipsJ','lh.slocs-v','rh.slocs-v']
    procSulci(['pmfs_p','aipsJ','slos1'],sulcus_plot_names)
