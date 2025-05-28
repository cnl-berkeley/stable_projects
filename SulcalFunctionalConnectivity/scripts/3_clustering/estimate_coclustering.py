import os,sys
from numpy import *
from matplotlib import *
import matplotlib.pyplot as plt
import pylab as P
import scipy.stats as stats


def checkSubModules(pipeline,tag,norm):
    pipel = 'orig'

    hemis = 'both'

    sub_ids = genfromtxt('/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-both/regs/namelist.csv'%(pipeline,pipel,tag,norm),dtype='str')
    n_subs = len(sub_ids)

    n_clust = []
    for s in arange(n_subs):
        savedir='/home/weiner/shakki/scripts-new/analysis/output/infomap/%s/%s/%s-%s/%s'%(pipeline,pipel,tag,norm,sub_ids[s])
        res_file = '%s/infomap-%s_%s_%s-fixed.txt'%(savedir,sub_ids[s],norm,hemis) # isolated nodes!
        res = genfromtxt(res_file,dtype='int',delimiter=',') # (n_roi,n_thrs)
        res = array(res,dtype='int')
      
        # smallest threshold with one cluster
        n_clusters = []
        for t in arange(shape(res)[1]):
            t_list = res[:,t].tolist()
            tmp = [x for x in t_list if t_list.count(x)>1]
            n_clusters.append(len(unique(tmp)))

        n_clust.append(max(n_clusters))
    print('min',min(n_clust))
    print('max',max(n_clust))
    print('mean',mean(n_clust))



def checkOne(pipeline,tag,norm):
    pipel = 'orig'

    savedir='/home/weiner/shakki/scripts-new/analysis/output/infomap/%s/%s/%s-%s/mean'%(pipeline,pipel,tag,norm)
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    hemis = 'both'

    sub_ids = genfromtxt('/home/weiner/shakki/scripts-new/analysis/output/setup/%s/%s/setup%s/%s-both/regs/namelist.csv'%(pipeline,pipel,tag,norm),dtype='str')
    n_subs = len(sub_ids)

    # define the srange of proportional thresholds to use for estimation
    thrs = arange(0.01,0.5,0.01) # as used in doInfomap.py!
    auc_thrs = thrs # arange(0.01,0.49,0.01) # corresponds to 0.05-0.5
    auc_min_ind = 0 #where(thrs > auc_thrs[0])[0][0]
    auc_max_ind0 = len(thrs) #where(thrs > auc_thrs[-1])[0][0]

    for s in arange(n_subs):
        savedir='/home/weiner/shakki/scripts-new/analysis/output/infomap/%s/%s/%s-%s/%s'%(pipeline,pipel,tag,norm,sub_ids[s])
        res_file = '%s/infomap-%s_%s_%s-fixed.txt'%(savedir,sub_ids[s],norm,hemis) # isolated nodes!
        res = genfromtxt(res_file,dtype='int',delimiter=',') # (n_roi,n_thrs)
        res = array(res,dtype='int')
        n_roi = shape(res)[0]
        n_thrs = shape(res)[1]
      
        # for each label, check min threshold it is clustered with at least one other label (and not zero)
        min_thrs = nan*ones((n_roi,1))
        for r in arange(n_roi):
            for t in arange(n_thrs):
                if (res[r,t] != 0) and len(where(res[:,t] == res[r,t])[0]) > 1:
                    min_thrs[r] = nanmin([min_thrs[r],t])
        min_thrs[isnan(min_thrs)] = auc_max_ind0

        # smallest threshold with one cluster
        one_cluster_ind = where(sum(res,0) == n_roi)[0] 
        if len(one_cluster_ind):      
            auc_max_ind = one_cluster_ind[0]
        else:
            auc_max_ind = auc_max_ind0 # 0.3

        # per label
        coi_auc = zeros((n_roi,n_roi))
        for i in arange(n_roi):
            auc_inds = arange(min_thrs[i],auc_max_ind,1).astype(int)
            n_thr = len(auc_inds)
            for t in arange(n_thr):
                for j in arange(n_roi):
                   if res[i,auc_inds[t]] == res[j,auc_inds[t]]:
                       coi_auc[i,j] = coi_auc[i,j] + 1
            if n_thr > 0:
                coi_auc[i,:] = 1.0*coi_auc[i,:]/n_thr
        res_file = '%s/infomap-%s_%s_%s-coi3-man.txt'%(savedir,sub_ids[s],norm,hemis)
        savetxt(res_file,coi_auc)

        # check if each label pair was clustered together at 5% of the thresholds
        coi_auc_bin = zeros((n_roi,n_roi))
        coi_auc_bin[coi_auc >= 0.05] = 1

        if s == 0:
            coi_all = zeros((n_subs,n_roi,n_roi))
        coi_all[s,:,:] = coi_auc_bin

    # across subjects, % each label pair clustered together
    coi_auc_perc = nanmean(coi_all,0)
    savedir='/home/weiner/shakki/scripts-new/analysis/output/infomap/%s/%s/%s-%s/mean'%(pipeline,pipel,tag,norm)
    res_file = '%s/infomap-%s_%s-coiperc4-man.txt'%(savedir,norm,hemis)
    savetxt(res_file,coi_auc_perc)

    calcAveCoclust(coi_auc_perc)

    stat(coi_all)


def calcAveCoclust(coi_auc_perc):
    npz = load('/home/weiner/shakki/scripts-new/visualization/matrix/plotting_order_man.npz')
    module_colors = npz['node_colors']
    labels = npz['labels']

    coi_auc_perc[eye(42) == 1] = nan
    
    tmp = module_colors[:,1]
    d = {ni: indi for indi, ni in enumerate(set(tmp))}
    numbers = [d[ni] for ni in tmp]

    noncluster_mask = zeros_like(coi_auc_perc)

    
    for n in unique(numbers):
        inds = [i for i in range(len(numbers)) if numbers[i] == n]
        for i in inds:
            for j in inds:
                 noncluster_mask[i,j] = 1

    print('mean', nanmean(coi_auc_perc))
    print('in-cluster mean', nanmean(coi_auc_perc[noncluster_mask == 1]))
    print('out-cluster mean', nanmean(coi_auc_perc[noncluster_mask == 0]))


def stat(coi_all):
    npz = load('/home/weiner/shakki/scripts-new/visualization/matrix/plotting_order_man.npz')
    module_colors = npz['node_colors']

    coi_all[:,eye(42) == 1] = nan
    
    tmp = module_colors[:,1]
    d = {ni: indi for indi, ni in enumerate(set(tmp))}
    numbers = [d[ni] for ni in tmp]

    noncluster_mask = zeros_like(squeeze(coi_all[0,:,:]))
    
    for n in unique(numbers):
        inds = [i for i in range(len(numbers)) if numbers[i] == n]
        for i in inds:
            for j in inds:
                 noncluster_mask[i,j] = 1

    n_subs = shape(coi_all)[0]
    arr = zeros((n_subs,2),dtype=float)
    for s in arange(n_subs):
        arr[s,0] = nanmean(coi_all[s,noncluster_mask == 1])
        arr[s,1] = nanmean(coi_all[s,noncluster_mask == 0])

    res = stats.ttest_rel(arr[:,0],arr[:,1],axis=0,alternative='two-sided')


if __name__ == '__main__':

    pipeline = '-bp'

    tags = ['-sareg-fdreg','-fdreg']
    for tag in tags:
        for norm in ['nonorm']:
            checkOne(pipeline,tag,norm)

