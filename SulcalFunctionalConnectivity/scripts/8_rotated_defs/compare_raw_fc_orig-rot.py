import os,sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def getOrigData(datadir,sub,label1,label2):
    fc_file = '%s/sub-%s/ses-%s/seed-fc_labels_hybrid.csv'%(datadir,sub,sub[-1])
    data = pd.read_csv(fc_file)
    column_names = np.array(data.columns)
    idx1 = np.argwhere(column_names == label1)
    idx2 = np.argwhere(column_names == label2)
    data = np.array(data)[idx1,idx2]

    return data


def getRotData(datadir,sub,label1,label2,perm):
    fc_file = '%s/rotated/sub-%s/ses-%s/seed-fc_labels_hybrid_perm%d.csv'%(datadir,sub,sub[-1],perm)
    data = pd.read_csv(fc_file)
    column_names = data.columns
    idx1 = np.argwhere(column_names == label1)
    idx2 = np.argwhere(column_names == label2)
    data = np.array(data)[idx1,idx2]

    return data


def doOneSulcus(datadir,subs,label1,label2,outdir,do_plot=False):

    # collect data
    n_reps = 1000
    all_orig = np.zeros((n_subs,1))
    all_rot = np.zeros((n_subs,n_reps))
    for s in np.arange(n_subs):
        try:
            all_orig[s] = getOrigData(datadir,subs[s],label1,label2)

            for p in np.arange(n_reps):
                all_rot[s,p] = getRotData(datadir,subs[s],label1,label2,p+1)
        except:
            print('error with %s'%subs[s])

    if do_plot:
        orig_vec = np.array(all_orig).squeeze()
        rot_vec = all_rot.flatten()
        arr = np.concatenate((rot_vec,orig_vec))
        df = pd.DataFrame(arr,columns=['value'])
        df['measure'] = np.array(['rot']*len(rot_vec) + ['orig']*len(orig_vec))
        sns.displot(df,x='value',hue='measure',stat='density',common_norm=False)
        plt.savefig('%s/dist_%s-%s.png'%(outdir,label1,label2))

    # average rotated
    rot_mean = np.nanmean(all_rot,axis=1)

    return all_orig,rot_mean


if __name__ == '__main__':

    outdir = '%s/density_plots'%os.getcwd()
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p','prts','lfms','aalf']
    lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
    labels = np.array(lpfc_labels + lpar_labels)

    datadir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fc_seed2seed-bp/relmatch'

    namelistfile = '/home/weiner/shakki/scripts-new/preprocessing/regress_motion/-bp/orig/fc-sareg-fdreg/namelist.csv'
    subs = np.loadtxt(namelistfile,dtype='str')
    n_subs = len(subs)
    print('n_subs',n_subs,subs)

    labels_sel = ['lh.painfs_any','lh.pmfs_a','lh.pmfs_i','lh.pmfs_p','lh.aipsJ','lh.slos1','rh.painfs_any','rh.pmfs_a','rh.pmfs_i','rh.pmfs_p','rh.aipsJ','rh.slos1']
    n_labels = len(labels_sel)
    arr_orig = np.zeros((n_labels,n_labels,n_subs))*np.nan
    arr_rot = np.zeros((n_labels,n_labels,n_subs))*np.nan
    for l1 in np.arange(n_labels):
        for l2 in np.arange(l1+1,n_labels,1):
            label1 = labels_sel[l1]
            label2 = labels_sel[l2]
            print(label1,label2)
            orig,rot = doOneSulcus(datadir,subs,label1,label2,outdir,False)
            arr_orig[l1,l2,:] = orig.squeeze()
            arr_rot[l1,l2,:] = rot.squeeze()

    arr_orig = arr_orig[np.isfinite(arr_orig)]
    arr_rot = arr_rot[np.isfinite(arr_rot)]

    # plot
    orig_vec = arr_orig.flatten()
    rot_vec = arr_rot.flatten()
    arr = np.concatenate((rot_vec,orig_vec))
    df = pd.DataFrame(arr,columns=['Correlation'])
    df['measure'] = np.array(['Random patches']*len(rot_vec) + ['Sulci']*len(orig_vec))
    g = sns.displot(df,x='Correlation',hue='measure',stat='density',common_norm=False)
    ax = plt.gca()
    ax.set_xlim(-1.25,1.25)
    plt.ylabel("Density") 
    g._legend.set_title('')
    plt.savefig('%s/dist_means-polished.png'%outdir)

    # means, ranges, sd
    print(np.nanmean(orig_vec),np.nanmin(orig_vec),np.nanmax(orig_vec),np.std(orig_vec))
    print(np.nanmean(rot_vec),np.nanmin(rot_vec),np.nanmax(rot_vec),np.std(rot_vec))
    
    # statistics
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.levene.html
    from scipy.stats import levene
    res = levene(orig_vec,rot_vec)
    print(res)
    #LeveneResult(statistic=857.695043783694, pvalue=1.0492364829886853e-175)
