import os,sys
from numpy import *
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import pylab as P
import bct
from sklearn import metrics
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix
import numpy as np
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import permutation_test_score


lpfc_labels = ['ifs','pimfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
lpar_labels = ['slocs_v','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
mpar_labels = []
plot_names = lpfc_labels + lpar_labels + mpar_labels

labels2 = []
for hemi in ['lh','rh']:
    for name in plot_names:
        labels2.append('%s.%s'%(hemi,name))


import bct

def norm(data_orig):
    data = zeros(shape(data_orig))
    print(shape(data))
    n_subs = len(data)
    for s in arange(n_subs):
         A = array(data_orig[s])
         data[s] = bct.weight_conversion(A,'normalize',copy=True) # scale to [0,1]

    return data


# https://stackoverflow.com/questions/41859613/how-to-obtain-reproducible-but-distinct-instances-of-groupkfold
def RandomGroupKFold_split(groups, n, seed=None):  # noqa: N802
    """
    Random analogous of sklearn.model_selection.GroupKFold.split.

    :return: list of (train, test) indices
    """
    groups = pd.Series(groups)
    ix = np.arange(len(groups))
    unique = np.unique(groups)
    np.random.RandomState(seed).shuffle(unique)
    result = []
    for split in np.array_split(unique, n):
        mask = groups.isin(split)
        train, test = ix[~mask], ix[mask]
        result.append((train, test))

    return result


def svm_pair(X,y):
    SVCClf = SVC(kernel = 'linear',gamma = 'scale', shrinking = False) #, decision_function_shape='ovo')
    
    scaler = MinMaxScaler()
    X = scaler.fit_transform(X.astype(float64))

    groups = tile(arange(43),2)

    # leave 3 subjects out cross-validation repeated 10 times
    split_reps = []
    for r in arange(10):
        split_reps.extend(RandomGroupKFold_split(groups, 40, seed=r))

    ytests = []
    ypreds = []
    for i, (train_ix, test_ix) in enumerate(split_reps):
        #print(f"Fold {i}:")
        #print(f"  Train: index={train_ix}")
        #print(f"  Test:  index={test_ix}")

        X_train,X_test = X[train_ix,:],X[test_ix,:]
        y_train,y_test = y[train_ix],y[test_ix]

        _ = SVCClf.fit(X_train,y_train)
        y_pred = SVCClf.predict(X_test)

        ytests += list(y_test)
        ypreds += list(y_pred)
        # Weights assigned to the features when kernel="linear".
        # These coefficients can be used directly as a crude type of feature importance score.
        # -- to get an idea of the relative importance of the features
        if i == 0:
            coefs = copy(SVCClf.coef_)
        else:
            coefs += SVCClf.coef_ # (n_classes, n_features)
    
    coefs = coefs/(i+1) # from sum to mean

#    acc = metrics.accuracy_score(ytests, ypreds)

    score_iris, perm_scores_iris, pvalue_iris = permutation_test_score(
        SVCClf, X, y, scoring="accuracy", cv=split_reps, n_permutations=100
    )

    conf_mat = confusion_matrix(ytests,ypreds)

    return conf_mat,coefs,score_iris,pvalue_iris

## sulcal discriminability
# pairwise excluding connectivity of the sulci to be classified (self-connections strong predictors)

npz = load('/home/weiner/shakki/scripts-new/analysis/output/setup/-bp/orig/setup-fdreg/nonorm-both/regs/data/allsub.npz')
#npz = load('/home/weiner/shakki/scripts-new/analysis/output/setup/-bp/orig/setup-sareg-fdreg/nonorm-both/regs/data/allsub.npz')
data = npz['data']
n_subs = shape(data)[0]
n_labels = shape(data)[1]

data = norm(data)

total_conf_mat = zeros((n_labels,n_labels))
total_f1_mat = zeros((n_labels,n_labels))
total_trial_counts = zeros((n_labels,n_labels))
total_coefs = nan
total_acc = zeros((n_labels,n_labels))
total_p = zeros((n_labels,n_labels))
for l1 in arange(len(labels2)):
    ind1 = labels2.index(labels2[l1])
    X1 = squeeze(data[:,ind1,:])
    y1 = [0]*n_subs

    for l2 in arange(l1+1,len(labels2)):
        print(l1, l2)
        ind2 = labels2.index(labels2[l2])
        X2 = squeeze(data[:,ind2,:])
        y2 = [1]*n_subs

        X = vstack((X1,X2))
        y = y1 + y2
        y = array(y)
        
        # exclude connections of the investigated sulci to avoid being driven by self-connections (set to zero)
        X_pruned = delete(X,[ind1,ind2],axis=1)

        conf_mat,coefs,acc,p = svm_pair(X_pruned,y)
        total_conf_mat[ind1,ind1] += conf_mat[0,0]
        total_conf_mat[ind1,ind2] += conf_mat[0,1]
        total_conf_mat[ind2,ind2] += conf_mat[1,1]
        total_conf_mat[ind2,ind1] += conf_mat[1,0]

        total_acc[ind1,ind2] = acc
        total_acc[ind2,ind1] = acc

        total_p[ind1,ind2] = p
        total_p[ind2,ind1] = p

        total_trial_counts[ind1,ind1] += 10
        total_trial_counts[ind1,ind2] += 10
        total_trial_counts[ind2,ind2] += 10
        total_trial_counts[ind2,ind1] += 10

        coefs_full = zeros((1,n_labels))
        inds = arange(n_labels)
        inds = delete(inds,[ind1,ind2],axis=0)
        coefs_full[0,inds] = coefs
        if (l1 == 0) and (l2 == 1):
            total_coefs = coefs_full
        else:
            total_coefs = vstack((total_coefs,coefs_full))

# scale by trials associated with each cell
total_conf_mat = divide(total_conf_mat,total_trial_counts)

savez('sulcal-classification-nosa-pair-scaled.npz',conf_mat=total_conf_mat,coefs=total_coefs,total_acc=total_acc,total_p=total_p)

