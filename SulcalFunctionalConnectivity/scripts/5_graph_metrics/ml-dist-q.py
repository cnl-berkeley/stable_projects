from numpy import *
import pandas as pd

# sa
npz = load('/home/weiner/shakki/scripts-new/analysis/group-analysis/#ml/sklearn/sulcal-classification-pair-scaled.npz')
acc_sa = npz['total_acc']

# nosa
npz = load('/home/weiner/shakki/scripts-new/analysis/group-analysis/#ml/sklearn/sulcal-classification-nosa-pair-scaled.npz')
acc_nosa = npz['total_acc']

# paired t test
from scipy.stats import ttest_rel
res = ttest_rel(acc_sa.flatten(),acc_nosa.flatten())
print(res,len(acc_sa)-1)
# -> accuracy overall higer with spatial autocorrelation fix (t = 3.81, p < 0.001) two-sided


# find pairs within < 5mm
files = '/home/weiner/shakki/scripts-new/analysis/output/setup/mean_adj.csv'
data = loadtxt(files,delimiter=',')
close_inds = where(data < 5)
acc_sa_close = acc_sa[close_inds].flatten()
acc_nosa_close = acc_nosa[close_inds].flatten()
res = ttest_rel(acc_sa_close,acc_nosa_close)
print(res,len(acc_sa_close)-1)
# -> accuracy higher with spatial autocorrelation fix for adjacent sulcal pairs (t = 2.43, p < 0.05) two-sided


# find crosshemispheric pairs
vals1 = []
vals2 = []
for i in arange(0,21,1):
    vals1.append(acc_sa[i,i+21])
    vals2.append(acc_nosa[i,i+21])
res = ttest_rel(vals1,vals2)
print(res,len(vals1)-1)
# -> accuracy higher with spatial autocorrelation fix also for crosshemispheric pairs (t = 2.98, p < 0.01)


# within hemisphere, > 10mm
close_inds = where(data > 10)
acc_sa[close_inds] = nan
acc_sa_lh = acc_sa[0:21,0:21].flatten()
acc_sa_rh = acc_sa[21:-1,21:-1].flatten()
acc_sa_whemi = concatenate((acc_sa_lh,acc_sa_rh))
acc_sa_whemi = acc_sa_whemi[~isnan(acc_sa_whemi)]

acc_nosa[close_inds] = nan
acc_nosa_lh = acc_nosa[0:21,0:21].flatten()
acc_nosa_rh = acc_nosa[21:-1,21:-1].flatten()
acc_nosa_whemi = concatenate((acc_nosa_lh,acc_nosa_rh))
acc_nosa_whemi = acc_nosa_whemi[~isnan(acc_nosa_whemi)]


res = ttest_rel(acc_sa_whemi,acc_nosa_whemi)
print(res,len(acc_sa_whemi)-1)
# -> accuracy not higher with spatial autocorrelation fix .5878834000701405, pvalue=0.11444681488808531)




lpfc_labels = ['ifs','pimfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
lpar_labels = ['slocs_v','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
mpar_labels = []
plot_names = lpfc_labels + lpar_labels + mpar_labels

labels2 = []
for hemi in ['lh','rh']:
    for name in plot_names:
        labels2.append('%s.%s'%(hemi,name))

inds = where(acc_sa < acc_nosa)
for i in arange(len(inds[0])):
   print(labels2[inds[0][i]],labels2[inds[1][i]])
