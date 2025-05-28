import os,sys
from numpy import *
import bct
import scipy.stats
from statsmodels.stats.multitest import multipletests
from copy import deepcopy

npz = load('/home/weiner/shakki/scripts-new/visualization/matrix/plotting_order_man.npz')
ci = npz['ciu']

thrs = arange(0.01,0.2,0.01)
n_thrs = len(thrs)

def compToMean(sub_vals):
    n_labels = shape(sub_vals)[1]
    sub_means = mean(sub_vals,axis=1)
    res = zeros((n_labels,2))
    for n in arange(n_labels):
        t,p = scipy.stats.ttest_ind(sub_vals[:,n],sub_means,alternative='two-sided')
        res[n,0] = p
        res[n,1] = t

    return res


def multiCorr(p_vals):
    res = multipletests(p_vals, alpha=0.05, method='fdr_bh') 
    return res[1]


def compAUCMeasure(data,savedir,metric_name,bin_tag):
    n_subs = shape(data)[0]
    n_labels = shape(data)[1]
    sub_vals = zeros((n_subs,n_labels,n_thrs))
    for s in arange(n_subs):
        for t in arange(n_thrs):
            data_thr = bct.threshold_proportional(data[s,:,:],thrs[t])
            sub_vals[s,:,t] = compMeasure(data_thr,metric_name,bin_tag)

    sub_vals2 = sum(sub_vals,2)
    savez('%s/auc-%s.npz'%(savedir,metric_name),data=sub_vals2)
    res = compToMean(sub_vals2)
    p = res[:,0]
    t = res[:,1]
    savetxt('%s/tvsmean-node_%s-auc.txt'%(savedir,metric_name),p)
    p_corr = multiCorr(squeeze(p).flatten())
    savetxt('%s/tvsmean-node_%s-auc_fdr.txt'%(savedir,metric_name),p_corr)
    savetxt('%s/tvsmean-node_t-%s-auc.txt'%(savedir,metric_name),t)

    mean_vals = mean(sum(sub_vals,2),0)
    savetxt('%s/mean-%s-auc.txt'%(savedir,metric_name),mean_vals)
    std_vals = std(sum(sub_vals,2),0)
    savetxt('%s/std-%s-auc.txt'%(savedir,metric_name),std_vals)


def compMeasure(data,metric_name,bin_tag):
    if bin_tag == 'bin':
        data = bct.weight_conversion(data,'binarize')
        if metric_name == 'pc':
            return bct.participation_coef(data,ci)
        elif metric_name == 'bw':
            return bct.betweenness_bin(data)
        elif metric_name == 'str':
            return bct.degrees_und(data)
    else:
        if metric_name == 'pc':
            return bct.participation_coef(data,ci)
        elif metric_name == 'bw':
            l = bct.weight_conversion(data,'lengths')
            return bct.betweenness_wei(l)
        elif metric_name == 'str':
            return bct.strengths_und(data)


if __name__ == '__main__':

    pipe = '-bp'

    pipeline = 'orig'

    bin_tag = 'bin'

    hemis = 'both'

    norm = 'nonorm'

    for tag in ['-fdreg','-sareg-fdreg']:
        subs = genfromtxt('/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-%s/regs/namelist.csv'%(pipe,pipeline,tag,norm,hemis),dtype='str')

        setupdir = '/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-%s/regs'%(pipe,pipeline,tag,norm,hemis)

        npz = load('%s/data/allsub.npz'%setupdir)
        data = npz['data']

        savedir='%s/graph-nopers-%s-tt'%(setupdir,bin_tag)
        print('saving to',savedir)
        if not os.path.exists(savedir):
            os.makedirs(savedir)

        ### whole network connectivity

        for metric_name in ['str','bw','pc']:
            compAUCMeasure(data,savedir,metric_name,bin_tag)
