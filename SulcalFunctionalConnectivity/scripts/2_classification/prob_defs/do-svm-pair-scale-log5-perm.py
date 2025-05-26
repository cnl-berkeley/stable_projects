import os,sys
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import pylab as P
import bct
from sklearn import metrics
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix
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

def norm(data_orig):
    data = np.nan*np.zeros(np.shape(data_orig))
    print(np.shape(data))
    n_subs = len(data)
    for s in np.arange(n_subs):
         A = np.array(data_orig[s])
         mask = ~np.isnan(A)
         A_norm = np.nan*np.zeros_like(A)
         A_norm[mask] = bct.weight_conversion(A[mask],'normalize',copy=True) # scale to [0,1]
         data[s] = A_norm

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
    X = scaler.fit_transform(X.astype(np.float64))

    groups = np.tile(np.arange(43),2)

    # leave 3 subjects out cross-validation repeated 10 times
    split_reps = []
    for r in np.arange(10):
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
            coefs = np.copy(SVCClf.coef_)
        else:
            coefs += SVCClf.coef_ # (n_classes, n_features)
    
    coefs = coefs/(i+1) # from sum to mean

#    acc = metrics.accuracy_score(ytests, ypreds)

    score_iris, perm_scores_iris, pvalue_iris = permutation_test_score(
        SVCClf, X, y, scoring="accuracy", cv=split_reps, n_permutations=100
    )
    
#    print('new acc',score_iris,pvalue_iris)

    conf_mat = confusion_matrix(ytests,ypreds)

    return conf_mat,coefs,score_iris,pvalue_iris


if __name__ == '__main__':
    npz_name = sys.argv[1]
    savename = sys.argv[2]
    outdir = '/home/weiner/shakki/scripts-new/analysis/group-analysis/ml-prob/sklearn_mpm'

    ## sulcal discriminability
    # pairwise excluding connectivity of the sulci to be classified (self-connections strong predictors)

    npz = np.load(npz_name)
    data = npz['data']
    n_subs = np.shape(data)[0]
    n_labels = np.shape(data)[1]

    # sulci whole connectivity values should not be included in the analysis
    excl_sulci = ['lh.slocs_v','lh.aipsJ','rh.pmfs_a','rh.aipsJ']
    excl_inds = []
    for excl_sulcus in excl_sulci:
        excl_ind = np.argwhere(np.array(labels2) == excl_sulcus)[0][0] # [10, 20, 23, 41]
        excl_inds.append(excl_ind)
        data[:,excl_ind,:] = np.nan
        data[:,:,excl_ind] = np.nan

    data = norm(data)

    total_conf_mat = np.zeros((n_labels,n_labels))
    total_f1_mat = np.zeros((n_labels,n_labels))
    total_trial_counts = np.zeros((n_labels,n_labels))
    total_coefs = np.nan
    total_acc = np.zeros((n_labels,n_labels))
    total_p = np.zeros((n_labels,n_labels))
    for l1 in np.arange(len(labels2)):
        ind1 = labels2.index(labels2[l1])
        X1 = np.squeeze(data[:,ind1,:])
        y1 = [0]*n_subs

        for l2 in np.arange(l1+1,len(labels2)):
            print(l1, l2)
            ind2 = labels2.index(labels2[l2])
            X2 = np.squeeze(data[:,ind2,:])
            y2 = [1]*n_subs

            X = np.vstack((X1,X2))
            y = y1 + y2
            y = np.array(y)
        
            # exclude connections of the investigated sulci to avoid being driven by self-connections (set to zero)
            X_pruned = np.delete(X,[ind1,ind2]+excl_inds,axis=1)

            try:
                conf_mat,coefs,acc,p = svm_pair(X_pruned,y)

                total_conf_mat[ind1,ind1] += conf_mat[0,0]
                total_conf_mat[ind1,ind2] += conf_mat[0,1]
                total_conf_mat[ind2,ind2] += conf_mat[1,1]
                total_conf_mat[ind2,ind1] += conf_mat[1,0]
                total_trial_counts[ind1,ind1] += 10
                total_trial_counts[ind1,ind2] += 10
                total_trial_counts[ind2,ind2] += 10
                total_trial_counts[ind2,ind1] += 10
            except:
                coefs = np.nan*np.zeros((1,X_pruned.shape[1]))
                p = np.nan
                acc = np.nan

            total_acc[ind1,ind2] = acc
            total_acc[ind2,ind1] = acc

            total_p[ind1,ind2] = p
            total_p[ind2,ind1] = p

            coefs_full = np.zeros((1,n_labels))
            inds = np.arange(n_labels)
            inds = np.delete(inds,[ind1,ind2]+excl_inds,axis=0)
            if len(inds) > 0:
                coefs_full[0,inds] = coefs
            if (l1 == 0) and (l2 == 1):
                total_coefs = coefs_full
            else:
                total_coefs = np.vstack((total_coefs,coefs_full))

    # scale by trials associated with each cell
    total_conf_mat = np.divide(total_conf_mat,total_trial_counts)
    print(total_conf_mat)

    # np.savez('%s/sulcal_classification_%s.npz'%(outdir,savename),conf_mat=total_conf_mat,coefs=total_coefs,total_acc=total_acc,total_p=total_p)

    total_acc[excl_inds,:] = np.nan
    total_acc[:,excl_inds] = np.nan
    np.savez('%s/sulcal_classification_%s_excl.npz'%(outdir,savename),conf_mat=total_conf_mat,coefs=total_coefs,total_acc=total_acc,total_p=total_p)


