import os,sys
from numpy import *
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import pylab as P
from statsmodels.stats.multitest import multipletests
import matplotlib as mpl

from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()

import cmasher as cmr

lpfc_labels = ['ifs','pimfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
mpar_labels = []
plot_names = lpfc_labels + lpar_labels + mpar_labels

labels = []
labels2 = []
for hemi in ['lh','rh']:
    for name in plot_names:
        labels.append('%s.%s'%(hemi,name))
        labels2.append('%s.%s'%(hemi,name))

n_names = len(labels)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(linspace(minval, maxval, n)))
    return new_cmap
cmap_mean = truncate_colormap(cmr.fusion_r,0.15,0.85)


def saveFig(vals_all,sulcus_plot_names,savedir,savename):
    fig = P.figure(figsize=(12,12))
    
    data = vals_all[:,:,0]

    a = fig.add_axes([0.25, 0.25, 0.05, 0.1])
    data[data == 0] = nan
    im = a.imshow(data, extent=[0,1,0,1], aspect='auto', interpolation='nearest', cmap=cmap_mean, origin='upper',vmin=-4,vmax=4,alpha=0.9)
    lims = im.get_clim()

    x = a.get_xlim()
    y = a.get_ylim()

    data_p = vals_all[:,:,1]
    data_p[isnan(data_p)] = 100
    print(shape(data_p))

    n_names2 = shape(vals_all)[0]
    xd = 1.0*(x[1]-x[0])/3
    yd = 1.0*(y[1]-y[0])/6
    for i in range(0,n_names2,1):
        for j in range(0,3,1):
            if i!=j and data_p[i][j] < 0.05:
            	a.plot(x[0]+0.5*xd+xd*j,y[1]-0.5*yd-yd*i,marker='o',markersize=4,color='k')

    # grid lines
    y_ticks = []
    for i in range(0,6,1): 
        y_ticks.append(y[1]-0.5*yd-yd*i)
    x_ticks = []
    for i in range(0,3,1):
        x_ticks.append(x[1]-0.5*xd-xd*i)

    a.spines['left'].set_linewidth(0.5)
    a.spines['right'].set_linewidth(0)
    a.spines['top'].set_linewidth(0)
    a.spines['bottom'].set_linewidth(0.5)

    xlabels = ['Degree','Betweenness','Participation']

    a.set_yticks(y_ticks)
    a.xaxis.tick_bottom()
    ylabels2 = sulcus_plot_names
    a.set_yticklabels(ylabels2,fontsize=10)
    a.set_xticks(x_ticks)
    a.set_xticklabels(xlabels[::-1],fontsize=10)
    a.xaxis.set_label_position('bottom')
    plt.xticks(rotation = 90)
    a.tick_params(axis='x', length=3)
    a.tick_params(axis='y', length=3)

    a = fig.add_axes([0.72,0.25,0.01,0.13])
    cbar = plt.colorbar(im,cax=a)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=0)
    cbar.set_ticks([lims[0],0,lims[1]])
    cbar.ax.set_yticklabels(['%.1f'%lims[0],'0','%.1f'%lims[1]],fontsize=10)
    cbar.set_label('$t$ value', labelpad=-10)

    if not os.path.exists(savedir):
        os.makedirs(savedir)

    fig.savefig('%s/%s.png'%(savedir,savename),bbox_inches='tight',pad_inches=0.03, dpi=300)


import statsmodels.api as sm

def regression(df,x_names,y_names):
    x = df[x_names]
    y = df[y_names]

    x = sm.add_constant(x)
    try:
        model = sm.OLS(y, x)
        res = model.fit()

        t = res.tvalues[x_names[0]]
        p = res.pvalues[x_names[0]]
    except:
        t = nan
        p = nan

    return t,p


def multiCorr(p_vals):
    res = multipletests(p_vals, alpha=0.05, method='fdr_bh') 
    return res[1]


def regMeasureNodal(measure,reg,tag,pipe,pipeline,bin_tag,savedir):

    if reg == 'MatrixRawScore':
        setupdir = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/nonorm-both/regs-MatrixRawScore'%(pipe,pipeline,tag)
    else:
        setupdir = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/nonorm-both/regs'%(pipe,pipeline,tag)

    infile = '%s/graph-nopers-%s/auc-%s.npz'%(setupdir,bin_tag,measure)
    npz = load(infile,allow_pickle=True)
    data = npz['data']

    df = pd.DataFrame(data,columns=labels2)
    
    df['sub'] = array(pd.read_csv('%s/namelist.csv'%setupdir,header=None))
    df[reg] = array(pd.read_csv('%s/%s.csv'%(setupdir,reg),header=None))
    df = df.melt(id_vars = ['sub',reg],value_vars = labels2,value_name = measure)

    vals_all = zeros((len(labels2),3))
    
    # correlations with the measure
    counter = 0
    for label in labels2:
        df_label = df[df.variable == label]
        vals_all[counter,0:2] = regression(df_label,[reg],[measure])
       
        counter = counter+1
    vals_all[:,2] = multiCorr(squeeze(vals_all[:,1]).flatten())

    savetxt('%s/reg-%s-%s_fdr.txt'%(savedir,measure,reg),squeeze(vals_all[:,2]))
    savetxt('%s/reg-%s-%s_uncp.txt'%(savedir,measure,reg),squeeze(vals_all[:,1]))
    savetxt('%s/reg-%s-%s_t.txt'%(savedir,measure,reg),squeeze(vals_all[:,0]))


def regMeasureGlobal(measure,reg,tag,pipe,pipeline,bin_tag,savedir):

    if reg == 'MatrixRawScore':
        setupdir = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/nonorm-both/regs-MatrixRawScore'%(pipe,pipeline,tag)
    else:
        setupdir = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/nonorm-both/regs'%(pipe,pipeline,tag)

    infile = '%s/graph-nopers-%s/auc-%s.npz'%(setupdir,bin_tag,measure)
    npz = load(infile,allow_pickle=True)
    data = npz['data']

    df = pd.DataFrame(data,columns=[measure])
    
    df['sub'] = array(pd.read_csv('%s/namelist.csv'%setupdir,header=None))

    df[reg] = array(pd.read_csv('%s/%s.csv'%(setupdir,reg),header=None)).flatten()

    vals_all = regression(df,[reg],[measure])
    print('efficiency %s'%reg,vals_all)


def regDepth(measure,sulcus,hemi,tag,pipe,pipeline,bin_tag,savedir):

    setupdir = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/nonorm-both/regs'%(pipe,pipeline,tag)

    infile = '%s/graph-nopers-%s/auc-%s.npz'%(setupdir,bin_tag,measure)
    npz = load(infile,allow_pickle=True)
    data = npz['data']

    df = pd.DataFrame(data,columns=labels2)
    
    df['sub'] = array(pd.read_csv('%s/namelist.csv'%setupdir,header=None))

    ind = labels2.index('%s.%s'%(hemi,sulcus))
    vals = squeeze(data[:,ind])
    reg = 'sulcal_depth_max_%s.%s'%(hemi,sulcus)
    df[reg] = array(pd.read_csv('%s/%s.csv'%(setupdir,reg),header=None)).flatten()

    df = df.melt(id_vars = ['sub',reg],value_vars = labels2,value_name = measure)

    # correlations with the measure
    df_label = df[df.variable == '%s.%s'%(hemi,sulcus)]
    vals_all = regression(df_label,[reg],[measure])
    print('%s sulcal_depth_max %s\t%s'%(measure,hemi,sulcus),vals_all)

    # plot scatter
    if measure == 'str':
        ylabel = 'Degree (%s.%s)'%(hemi,sulcus)
    else:
        ylabel = '%s (%s.%s)'%(measure,hemi,sulcus)

    xlabel = 'Depth (%s.%s)'%(hemi,sulcus)
    plotScatter(df_label,measure,reg,'%s-%s'%(measure,reg),xlabel,ylabel,savedir)

    return vals_all


def plotScatter(df_label,col1,col2,savename,ylabel,xlabel,savedir):

    fig = P.figure(figsize=(5,5))

    a = fig.add_axes([0.25, 0.25, 0.46, 0.46])

    sns.regplot(x=col1,y=col2,data=df_label,ax=a, scatter_kws={'s':15})

    a.set_xlabel(xlabel, size = 9)
    a.set_ylabel(ylabel, size = 9)
    a.spines['left'].set_linewidth(0.5)
    a.spines['right'].set_linewidth(0)
    a.spines['top'].set_linewidth(0)
    a.spines['bottom'].set_linewidth(0.5)

    plt.savefig('%s/scatter-%s.png'%(savedir,savename),bbox_inches='tight')
    plt.clf()


def regDepthWAge(measure,sulcus,hemi,tag,pipe,pipeline,bin_tag,savedir):

    setupdir = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/nonorm-both/regs'%(pipe,pipeline,tag)

    infile = '%s/graph-nopers-%s/auc-%s.npz'%(setupdir,bin_tag,measure)
    npz = load(infile,allow_pickle=True)
    data = npz['data']

    df = pd.DataFrame(data,columns=labels2)
    
    df['sub'] = array(pd.read_csv('%s/namelist.csv'%setupdir,header=None))

    df['age'] = array(pd.read_csv('%s/ScanAge.csv'%setupdir,header=None))

    ind = labels2.index('%s.%s'%(hemi,sulcus))
    vals = squeeze(data[:,ind])
    reg = 'sulcal_depth_max_%s.%s'%(hemi,sulcus)
    df[reg] = array(pd.read_csv('%s/%s.csv'%(setupdir,reg),header=None)).flatten()
    #df['mFD'] = array(pd.read_csv('%s/mFD.csv'%setupdir,header=None))

    df['age'] = df['age'].astype(float)
    df[reg] = df[reg].astype(float)
    df['age'] = scaler.fit_transform(df['age'].values.reshape(-1,1))
    df[reg] = scaler.fit_transform(df[reg].values.reshape(-1,1))

    df = df.melt(id_vars = ['sub','age',reg],value_vars = labels2,value_name = measure)

    # correlations with the measure
    df_label = df[df.variable == '%s.%s'%(hemi,sulcus)]
    vals_all = regression(df_label,[reg,'age'],[measure])
    print('%s sulcal_depth_max wage %s\t%s'%(measure,hemi,sulcus),vals_all)

    return vals_all


def procSulci(sulci,sulcus_plot_names):

    pipe = '-bp'

    pipeline = 'orig'

    tag = '-sareg-fdreg' #'-fdreg'

    bin_tag = 'bin'

    savedir = '/home/weiner/shakki/scripts-new/analysis/output/visualization/node_regs/%s/%s/%s/%s'%(bin_tag,pipe,tag,pipeline)
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    measures = ['str','bw','pc']
    res_all = zeros((len(sulci)*2,len(measures),4))*nan # sulcus, measure, (t,p,p_fdr_sulcus,p_fdr_all)
    counter = 0
    for sulcus in sulci:
        for hemi in ['lh','rh']:
            for m in arange(len(measures)):
                res_all[counter,m,0:2] = regDepthWAge(measures[m],sulcus,hemi,tag,pipe,pipeline,bin_tag,savedir)
            res_all[counter,:,2] = multiCorr(squeeze(res_all[counter,:,1]))
            counter = counter+1
    p_vals = squeeze(res_all[:,:,1]).flatten()
    print(shape(p_vals))
    res_all[:,:,3] = reshape(multiCorr(p_vals),(len(sulci)*2,len(measures)))

    set_printoptions(formatter={'float': '{: 0.3f}'.format})

    print(res_all)

    res_all[:,:,3] = reshape(multiCorr(p_vals),(6,3))

    # make colormesh plot
    res_all = res_all[:,:,[0,3]]
    print(shape(res_all))
    saveFig(res_all,sulcus_plot_names,savedir,'metrics-depth-mesh_%s'%sulci[0])


if __name__ == '__main__':
    sulcus_plot_names = ['lh.pimfs','rh.pimfs','lh.pmfs-a','rh.pmfs-a','lh.pmfs-i','rh.pmfs-i'] # only for plotting
    procSulci(['pimfs_any','pmfs_a','pmfs_i'],sulcus_plot_names)

    sulcus_plot_names = ['lh.pmfs-p','rh.pmfs-p','lh.aipsJ','rh.aipsJ','lh.slocs-v','rh.slocs-v']
    procSulci(['pmfs_p','aipsJ','slos1'],sulcus_plot_names)
