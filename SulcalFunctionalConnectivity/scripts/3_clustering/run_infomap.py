import os,sys
from numpy import *
from matplotlib import *
import matplotlib.pyplot as plt
import pylab as P
import networkx as nx
from infomap import Infomap
import bct

def getLabels(hemis):
    lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
    lpar_labels = ['slocs1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    mpar_labels = []
    plot_names = lpfc_labels + lpar_labels + mpar_labels

    labels = []
    for hemi in hemis:
        for name in plot_names:
            labels.append('%s.%s'%(hemi,name))

    return labels


def do(A,savedir,savename,labels):
    A[A < 0] = 0
    print('* shape(A)',shape(A))

    G = nx.from_numpy_array(A)
    n_nodes = len(labels)

    # define thresholds to check
    thrs = arange(0.01,0.5,0.01)

    # calculate clustering at each threshold
    res = zeros((n_nodes,len(thrs)))
    t_inds = arange(len(thrs)).tolist()
    for t in arange(len(t_inds)):
        A_thr = bct.threshold_proportional(A,thrs[t_inds[t]],copy=True)
        G_thr = nx.from_numpy_array(A_thr)
        G_thr.edges(data=True)
        im = Infomap(silent=False,seed=123)
        mapping = im.add_networkx_graph(G_thr,weight='weight')
        im.run()

        for node in im.nodes:
            res[node.node_id,t_inds[t]] = node.module_id

    res_path = '%s/infomap-%s.txt'%(savedir,savename)
    savetxt(res_path,res)
    plot(res,thrs,t_inds[9],'%s/infomap-%s.png'%(savedir,savename))

    # try to sort colors so clusters would match across thresholds
    cmd="matlab-2021a -nodisplay -nosplash -nodesktop -r \"ciu_value_fix('%s'); exit\""%res_path
    os.system(cmd)
    res_fixed_path = '%s/infomap-%s-fixed.txt'%(savedir,savename)
    res_fixed = genfromtxt(res_fixed_path,dtype='int',delimiter=',')
    plot(res_fixed,thrs,t_inds[9],'%s/infomap-%s-fixed.png'%(savedir,savename))


def plot(res,thrs,thr_for_sorting,savepath):
    # order data for visualization
    order = argsort(res[:,thr_for_sorting])
    res_ordered = res[order,:]
    y_ordered = array(labels)[order]

    # plot
    plt.rcParams['axes.facecolor']='white'
    plt.rcParams['savefig.facecolor']='white'
    fig = P.figure(figsize=(10,5))
    a = fig.add_axes([0.05, 0.05, 0.9, 0.9])
    im = a.imshow(res_ordered,cmap = plt.cm.jet)
    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.set_xticks(arange(len(thrs)))
    n_nodes = shape(res)[0]
    a.set_yticks(arange(n_nodes))
    a.set_xticklabels(thrs)
    plt.xticks(rotation = 90)
    a.set_yticklabels(y_ordered)

    fig.savefig(savepath,bbox_inches='tight',pad_inches=0.03, dpi=300)
    plt.close()


if __name__ == '__main__':

    pipe = '-bp'

    pipeline = 'orig'

    tags = ['-fdreg','-sareg-fdreg']
    for tag in tags:

        # mean networks
        for norm in ['nonorm']:
            savedir='/home/weiner/shakki/scripts-new/analysis/output/infomap/%s/%s/%s-%s/mean'%(pipe,pipeline,tag,norm)
            if not os.path.exists(savedir):
                os.makedirs(savedir)

            for hemis in ['both']: #,'lh','rh']:
                datafile = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-%s/regs/data/mean.csv'%(pipe,pipeline,tag,norm,hemis)
                A = loadtxt(datafile,delimiter=',')

                if hemis == 'both':
                    labels = getLabels(['lh','rh'])
                else:
                    labels = getLabels([hemis])
                savename = 'mean_%s_%s'%(norm,hemis)
                do(A,savedir,savename,labels)

        # networks of individual subjects
        npz = load('/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-%s/regs/data/allsub.npz'%(pipe,pipeline,tag,norm,hemis))
        sub_ids = genfromtxt('/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-%s/regs/namelist.csv'%(pipe,pipeline,tag,norm,hemis),dtype='str')
        data = npz['data']
        n_subs = shape(data)[0]

        for hemis in ['both']:
            if hemis == 'both':
                labels = getLabels(['lh','rh'])
            else:
                labels = getLabels([hemis])

            for norm in ['nonorm','posnorm']:
                for s in arange(n_subs):
                    savedir='/home/weiner/shakki/scripts-new/analysis/output/infomap/%s/%s/%s-%s/%s'%(pipe,pipeline,tag,norm,sub_ids[s])
                    if not os.path.exists(savedir):
                        os.makedirs(savedir)
                    A = squeeze(data[s,:,:])
                    savename = '%s_%s_%s'%(sub_ids[s],norm,hemis)
                    do(A,savedir,savename,labels)

