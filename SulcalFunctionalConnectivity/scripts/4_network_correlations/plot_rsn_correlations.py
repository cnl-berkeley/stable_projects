import seaborn as sns
import pandas as pd
import pylab as P
import matplotlib.pyplot as plt
import os,sys
from numpy import *
from matplotlib import *
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
import numpy as np
import cmasher as cmr

npz = load('out/plotting_order_man_pTS.npz')
order = npz['order']
module_colors = npz['node_colors']
labels = npz['labels']
c = npz['c']

n_names = len(labels)
labels_ordered = array(labels)[order]
labels_ordered[labels_ordered == 'lh.painfs_any'] = 'lh.pimfs'
labels_ordered[labels_ordered == 'rh.painfs_any'] = 'rh.pimfs'

def getLabelInds():
    # wanted label names
    lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
    lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    mpar_labels = []    
    req_labels = lpfc_labels + lpar_labels + mpar_labels

    # label names supplied to the script to calculate seed2annot connectivity
    lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
    lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    mpar_labels = []

    all_labels = lpfc_labels + lpar_labels + mpar_labels

    # corresponding indices in the seed2annot output
    label_inds = []
    for i in arange(len(all_labels)):
        if all_labels[i] in req_labels:
            label_inds.append(i)

    return label_inds,req_labels


def pruneFC(datadir,sub,hemi,label_inds):
    # select data for wanted labels and networks

    fc_file = '%s/sub-%s/ses-%s/hemi-%s/seed-fc-annot_%s_nw_hybrid.csv'%(datadir,sub,sub[-1],hemi,hemi)
    print(fc_file)
    data = pd.read_csv(fc_file)
    print(data)
    column_names = data.columns
    print(column_names)
    data = array(data)
    print(shape(data))
    fc_arr = data[label_inds,1::]
    column_names = column_names[1::].tolist()
    for c in arange(len(column_names)):
       column_names[c] = column_names[c][2:-1]
    
    return fc_arr,column_names


def doOne(data,savepath,thrs,nw_names,cb_pars=None,mymap=cmr.fusion):
    plt.rcParams.update({'mathtext.default':'regular' })

    fig = P.figure(figsize=(12,12))
    fig,im = plotScalar(fig,data,thrs,0,0,nw_names,cb_pars,mymap)
 
    fig.savefig('%s.svg'%savepath, format="svg",bbox_inches='tight',pad_inches=0.03, dpi=600)
    plt.close()


def plotScalar(fig,data_both,lims,row_ind,col_ind,nw_names,cb_pars,mymap):
    font0 = FontProperties(size=4)
    mx=0.2

    # sort x
    data = np.array(data_both[0])
    data = data[order]
    data = data.T

    datap = np.array(data_both[1])
    datap = datap[order]
    datap = datap.T

    # sort y
    h = 0.46/42*14
    a = fig.add_axes([0.25, 0.25, 0.46, h])

    data[data == 0] = nan
    lims = cb_pars[0]
    im = a.imshow(data, extent=[0,1,0,1], aspect='auto', origin='upper', interpolation='nearest', cmap=mymap, norm=mpl.colors.TwoSlopeNorm(lims[1],vmin=lims[0],vmax=lims[2]))

    x = a.get_xlim()
    y = a.get_ylim()

    xd = 1.0*(x[1]-x[0])/n_names
    yd = 1.0*(y[1]-y[0])/len(nw_names)
    
    # grid lines
    x_ticks = []
    for i in range(0,n_names,1):
        x_ticks.append(x[1]-0.5*xd-xd*i)
    y_ticks = []
    for i in range(0,len(nw_names),1):
        y_ticks.append(y[1]-0.5*yd-yd*i)

    # significance markings
    for i in range(0,len(nw_names),1):
        for j in range(0,n_names,1):
            if datap[i,j] < 0.05/2:
           	    a.plot(x[0]+0.5*xd+xd*j,y[1]-0.5*yd-yd*i,marker='o',markersize=4,color='k')

    a.spines['left'].set_linewidth(0.5)
    a.spines['right'].set_linewidth(0.0)
    a.spines['top'].set_linewidth(0.0)
    a.spines['bottom'].set_linewidth(0.5)

    a.set_yticks(y_ticks)
    a.xaxis.tick_bottom()
    a.set_yticklabels(nw_names)
    a.set_xticks(x_ticks)
    a.set_xticklabels(labels_ordered[::-1])
    a.xaxis.set_label_position('top')
    plt.xticks(rotation = 90)
    a.tick_params(axis='x', length=13)
    a.tick_params(axis='y', length=3)
    
    a = fig.add_axes([0.72,0.25,0.008,0.1])
    cbar = plt.colorbar(im,cax=a)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=0)
    if cb_pars == None:
        cbar.set_ticks([lims[0],0,lims[1]])
        cbar.ax.set_yticklabels(['%.2f'%lims[0],'0','%.2f'%lims[1]]) 
    else:
        cbar.set_ticks(cb_pars[0])
        cbar.ax.set_yticklabels(cb_pars[1])   
        cbar.set_label(cb_pars[2])

    a = fig.add_axes([0.25+0.0005, 0.25-0.0125, 0.46, 0.01]) # x orientation
    a.imshow(swapaxes(c,0,1), interpolation='nearest', aspect='auto')
    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks([])
    a.set_yticks([])

    return fig,im


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(linspace(minval, maxval, n)))
    return new_cmap


if __name__ == '__main__':

    datadir = '/home/weiner/shakki/scripts-github/analysis/network_correlations'

    label_inds,labels = getLabelInds()
    n_labels = len(labels)
    print('n_labels',len(label_inds))

    nw_names_full = {'FP':'Frontal parietal','DAN':'Dorsal attention','DMN':'Default mode','Sal':'Salience','CO':'Cingulo-opercular','VAN':'Ventral attention','PMN':'Parietal medial','Vis':'Visual','SMd':'Sensorimotor dorsal','PON':'Parieto-occipital','Aud':'Auditory','SMl':'Sensorimotor lateral','MTL':'Medial temporal','Tpole':'Temporal pole'}

    hemi = 'lh'
    data_all = np.zeros((n_labels*2,14))*np.nan
    p_all = np.zeros((n_labels*2,14))*np.nan
    for i in np.arange(n_labels):
        npz = np.load('%s/bootstrap/collected_results_%s.%s.npz'%(datadir,hemi,labels[i]),allow_pickle=True)
        collected_results = npz['arr_0']
        data_all[i,:] = collected_results[:,0] 
        p_all[i,:] = np.minimum(collected_results[:,2],collected_results[:,3])
        nw_names = npz['arr_1']
    hemi = 'rh'
    for i in np.arange(n_labels):
        npz = np.load('%s/bootstrap/collected_results_%s.%s.npz'%(datadir,hemi,labels[i]),allow_pickle=True)
        collected_results = npz['arr_0']
        data_all[n_labels+i,:] = collected_results[:,0]
        p_all[n_labels+i,:] = np.minimum(collected_results[:,2],collected_results[:,3])

    for key in nw_names_full:
        nw_names[nw_names == key] = nw_names_full[key]

    # sort y
    tmp = sum(data_all,axis=0)
    print(shape(tmp))
    yorder = np.argsort(tmp)[::-1]
    data_mean = data_all[:,yorder]
    data_p = p_all[:,yorder]
    nw_names = np.array(nw_names)
    nw_names = nw_names[yorder]

    cmap_var = cmr.lavender
    cmap_mean = truncate_colormap(cmr.fusion_r,0.15,0.85)
    
    outdir = '%s/out/matrix'%os.getcwd()
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    savepath = '%s/network_correlations'%outdir
    cb_pars = [[0,10,20],['0','10','20'],'t']
    print(data_mean)
    print(np.amax(data_mean),np.amin(data_mean))
    doOne([data_mean,data_p],savepath,[0,1],nw_names,cb_pars,cmap_mean)
