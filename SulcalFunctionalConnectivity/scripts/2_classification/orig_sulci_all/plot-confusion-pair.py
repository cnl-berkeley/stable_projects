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

from statsmodels.stats.multitest import multipletests

npz = load('/home/weiner/shakki/scripts-new/visualization/matrix/plotting_order_man.npz')
order = npz['order']
module_colors = npz['node_colors']
labels = npz['labels']
c = npz['c']
n_names = len(labels)
labels_ordered = array(labels)[order]

tert_sulci = ['pmfs-p','pmfs-i','pmfs-a','pimfs','aalf','prts','lfms','aipsJ','slocs-v']
label_type = []
for label in labels:
    if label.split('.')[1] in tert_sulci:
        label_type.append('tertiary')
    else:
        label_type.append('primary')


def plotScalar(fig,data,lims,row_ind,col_ind,mymap,cb_pars):
    font0 = FontProperties(size=4)
    mx=0.2

    data[data < lims[0]] = lims[0]+0.0001
    data[data > lims[1]] = lims[1]-0.0001

    data = array(data).astype(float)
    data = data[order][:,order]

    a = fig.add_axes([0.25, 0.25, 0.46, 0.46])

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
        cbar.set_label(cb_pars[2], labelpad=-15)

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


def plotScalarSig(fig,data,lims,row_ind,col_ind,mymap,cb_pars):
    font0 = FontProperties(size=4)
    mx=0.2

    data[0][data[0] < lims[0]] = lims[0]+0.0001
    data[0][data[0] > lims[1]] = lims[1]-0.0001

    data_v = array(data[0]).astype(float)
    data_v = data_v[order][:,order]

    a = fig.add_axes([0.25, 0.25, 0.46, 0.46])

    data_v = tril(data_v)
    data_v[data_v == 0] = nan
    im = a.imshow(data_v, origin='upper',extent=[0,1,0,1], aspect='auto', interpolation='nearest', cmap=mymap, vmin=lims[0], vmax=lims[1])

    x = a.get_xlim()
    y = a.get_ylim()

    xd = 1.0*(x[1]-x[0])/n_names
    yd = 1.0*(y[1]-y[0])/n_names
    
    # grid lines
    linemarks = arange(0.5,n_names+0.5,1)
    y_ticks = []
    for i in range(0,n_names,1):
        y_ticks.append(y[1]-0.5*yd-yd*i)

    ### p values
    p_vals = tril(data[1]).flatten()
    p_vals = p_vals[p_vals > 0]
    res = multipletests(p_vals, alpha=0.05, method='fdr_bh')[0]
    counter = 0
    data_p_fdr = zeros_like(data[1])
    for i in arange(n_names):
        for j in arange(i+1,n_names,1):
            if res[counter] == True:
                data_p_fdr[i,j] = 100
            counter = counter+1
    data_p_fdr = maximum(data_p_fdr,data_p_fdr.T)
    data_p_fdr = data_p_fdr[order][:,order]
    data_p_fdr = tril(data_p_fdr)

    print('check')
    print(unique(data_p_fdr))
    
    xd = 1.0*(x[1]-x[0])/n_names
    yd = 1.0*(y[1]-y[0])/n_names
    
    for i in range(0,n_names,1):
        for j in range(0,i,1):
            if i!=j and data_p_fdr[i][j] > 0:
            	a.plot(x[0]+0.5*xd+xd*j,y[1]-0.5*yd-yd*i,marker='o',markersize=4,color='r')

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
        cbar.set_label(cb_pars[2], labelpad=-10)

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


def plotSel(fig,data,ylabels,lims,row_ind,col_ind,mymap,cb_pars,plot_labels=True):
    font0 = FontProperties(size=4)
    mx=0.2

    data[data < lims[0]] = lims[0]+0.0001
    data[data > lims[1]] = lims[1]-0.0001

    for i in arange(len(ylabels)):
        ylabels[i] = ylabels[i].replace('lh.','')
        ylabels[i] = ylabels[i].replace('rh.','')

    data = array(data).astype(float)
    data = data[:,order]

    a = fig.add_axes([0.25, 0.25, 0.46, 0.46*3/42])

    im = a.imshow(data, origin='upper',extent=[0,1,0,1], aspect='auto', interpolation='nearest', cmap=mymap, vmin=lims[0], vmax=lims[1])

    x = a.get_xlim()
    y = a.get_ylim()

    xd = 1.0*(x[1]-x[0])/n_names
    yd = 1.0*(y[1]-y[0])/len(ylabels)
    
    # grid lines
    linemarks = arange(0.5,n_names+0.5,1)
    x_ticks = []
    for i in range(0,n_names,1):
        x_ticks.append(x[1]-0.5*xd-xd*i)

    linemarks = arange(0.5,len(ylabels)+0.5,1)
    y_ticks = []
    for i in range(0,len(ylabels),1):
        y_ticks.append(y[1]-0.5*yd-yd*i)

    a.spines['left'].set_linewidth(0.5)
    a.spines['right'].set_linewidth(0)
    a.spines['top'].set_linewidth(0)
    a.spines['bottom'].set_linewidth(0.5)

    a.set_yticks(y_ticks)
    a.xaxis.tick_bottom()
    a.set_xticks(x_ticks)
    a.set_yticklabels(ylabels)
    a.xaxis.set_label_position('bottom')
    plt.xticks(rotation = 90)
    a.tick_params(axis='y', length=3)
    if plot_labels:
        a.set_xticklabels(labels_ordered[::-1])
        a.tick_params(axis='x', length=13)

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
            cbar.set_label(cb_pars[2], labelpad=-20)
    else:
        a.set_xticklabels([])

    if plot_labels:
        a = fig.add_axes([0.25+0.0005, 0.25-0.0125, 0.46, 0.01]) # x orientation
        a.imshow(swapaxes(c,0,1), interpolation='nearest', aspect='auto')
        a.spines['left'].set_color('none')
        a.spines['right'].set_color('none')
        a.spines['top'].set_color('none')
        a.spines['bottom'].set_color('none')
        a.set_xticks([])
        a.set_yticks([])

    return fig,im


def doOne(data,savepath,thrs,cmap,ylabels=None,cb_pars=None):
    plt.rcParams['axes.facecolor']='white'
    plt.rcParams['savefig.facecolor']='white'
    plt.rcParams.update({'mathtext.default':'regular' })

    fig = P.figure(figsize=(12,12))
    if ylabels == None and len(data) == 2:
        fig,im = plotScalarSig(fig,data,thrs,0,0,cmap,cb_pars)
    elif ylabels == None:
        fig,im = plotScalar(fig,data,thrs,0,0,cmap,cb_pars)
    else:
        fig,im = plotSel(fig,data,ylabels,thrs,0,0,cmap,cb_pars)
    print('* saving')
    fig.savefig('%s.png'%savepath,bbox_inches='tight',pad_inches=0.03, dpi=300)
    plt.close()


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(linspace(minval, maxval, n)))
    return new_cmap


def plotTable(num_data,savepath):
    plt.rcParams['axes.facecolor']='white'
    plt.rcParams['savefig.facecolor']='white'
    plt.rcParams.update({'mathtext.default':'regular' })

    fig = P.figure(figsize=(12,12))
    a = fig.add_axes([0.25, 0.25, 0.46, 0.46])    
    
    num_data = num_data[:,order]
    num_data = around(num_data,decimals=3)
    print(shape(num_data))

    str_data = num_data.astype('str')
    cell_text = vstack((labels_ordered,str_data)).T
    columns = ['','Precision','Recall','F1 score']
    t = a.table(cellText=cell_text,colLabels=columns,colWidths=[.2,.2,.2,.2],loc='left',edges='open') #cellColours=colors,
    t.set_fontsize(16)
    for r in range(0, len(columns)):
        cell = t[0, r]
        cell.set_height(0.05)

    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks([])
    a.set_yticks([])

    print('* saving')
    fig.savefig('%s.png'%savepath,bbox_inches='tight',pad_inches=0.03, dpi=300)
    plt.close()


def getFeatureInds(slabel,tlabels,task):
    ylabels = []
    task_inds = []
    mult = [1]*len(tlabels)
    for t in arange(len(tlabels)):
        if tlabels[t] == slabel:
            task_inds.append(nan)
        else:
            try:
                tmp1 = where(task == '%s-%s'%(slabel,tlabels[t]))[0][0]
            except:
                tmp1 = -1
            try:
                tmp2 = where(task == '%s-%s'%(tlabels[t],slabel))[0][0]
                mult[t] = -1
            except:
                tmp2 = -1
            task_inds.append(max(tmp1,tmp2))
            
        ylabels.append('%s vs %s'%(slabel,tlabels[t]))
    task_inds = array(task_inds).flatten().astype(int)

    return task_inds,ylabels,mult


if __name__ == '__main__':

#    npz = load('sulcal-classification-pair-scaled.npz',allow_pickle=True)
    npz = load('sulcal-classification-nosa-pair-scaled.npz',allow_pickle=True)
    conf_mat = npz['conf_mat']

    total_p = npz['total_p']

    savedir=os.getcwd() #'/home/weiner/shakki/scripts-new/analysis/output/visualization/matrix/%s/%s'%(pipe,pipeline)
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    
    # plot confusion matrix
    cmap_var = truncate_colormap(cmr.rainforest_r,0,0.95)
    savepath = '%s/conf'%savedir
    doOne(conf_mat,savepath,[0,43],cmap_var,None,[[0,43],['0','100'],'Percent'])

    # plot accuracy matrix (from binary classifications)
    total_acc = npz['total_acc']
    total_acc[eye(shape(total_acc)[0]) == 1] = nan
    savepath = '%s/acc'%savedir
    doOne(total_acc,savepath,[0.5,1.0],cmr.rainforest,None,[[0.5,1.0],['50','100'],'Accuracy (%)'])
    # with permuted significance
    total_acc[eye(shape(total_acc)[0]) == 1] = nan
    #cmap_mean = truncate_colormap(cmr.fusion_r,0.15,0.85)
    savepath = '%s/accsig'%savedir
    doOne([total_acc,total_p],savepath,[0.5,1.0],cmr.rainforest,None,[[0.5,1.0],['50','100'],'Accuracy (%)'])

    print('Mean: %f (min %f, max %f)'%(nanmean(total_acc), nanmin(total_acc),nanmax(total_acc)))
