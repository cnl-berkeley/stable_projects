import os,sys
from numpy import *

def find_nearest_ind(array,value):
    ind = min(range(len(array)), key=lambda i: abs(array[i]-value))

    return ind 


def getWeightMatrix(adj,var_x,var_y):
    w_arr = zeros(shape(adj))
    n_roi = shape(adj)[0]
    for i in arange(n_roi):
       for j in arange(n_roi):
            ni = find_nearest_ind(var_x,adj[i,j])
            w_arr[i,j] = var_y[ni]

    return w_arr


def getResiduals(data,w_arr):
    n_roi = shape(data)[0]

    y = data.flatten()
    x = w_arr.flatten()

    # exclude self-connections
    i_ex = argwhere(eye(n_roi).flatten() == 1).flatten().tolist()
    y = delete(y,i_ex)
    x = delete(x,i_ex)

    # we are only interested in reducing the positive bias in short-distance connections,
    # so rescale the regressor so it will not modify longer distance connections
    w_arr = w_arr-amin(w_arr)

    # exclude nan connections (in control analysis with probabilistic seed labels)
    inds_defined = argwhere((y > 0) & (y != inf)).flatten()
    y = y[inds_defined]
    x = x[inds_defined]

    # solve the regression model
    from scipy import stats
    ret = stats.linregress(x,y)

    # collect residualized values per connection
    res_arr = zeros((n_roi,n_roi))
    for r1 in arange(n_roi):
       for r2 in arange(n_roi):
           print(r1,r2,data[r1,r2],w_arr[r1,r2],ret[0])
           res_arr[r1,r2] = data[r1,r2] - ret[0]*w_arr[r1,r2] # leave intercept

    return res_arr


if __name__ == '__main__':

    sub = sys.argv[1]
    hemi = sys.argv[2]
    pipeline = sys.argv[3]
    prob_tag = sys.argv[4]

    datadir='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation'

    # get distance matrix
    if prob_tag != 'prob':
        datadir_func = '/home/weiner/shakki/NotBackedUp/spatial_autocorrelation%s/%s/labelnw'%(pipeline,sub)
        adj_file = '/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/labelnw/adj-labels-%s.txt'%(sub,hemi)
    else:
        datadir_func = '/home/weiner/shakki/NotBackedUp/spatial_autocorrelation%s/%s/labelnw_prob_mpm'%(pipeline,sub)
        adj_file = '/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/labelnw_prob_mpm/adj-labels-%s.txt'%(sub,hemi)
    adj = loadtxt(adj_file)
    adj = adj*10 # convert to same units as fit_x; was in units of 1mm
    print(shape(adj))

    # get estimated autocorrelation
    try:
        var_x = loadtxt('%s/%s/cortex/fit_x-%s.txt'%(datadir,sub,hemi)) # max 200 mm i.e. 20 cm -> units of 10mm
        var_y = loadtxt('%s/%s/cortex/fit_y-%s.txt'%(datadir,sub,hemi))
    except: # use the estimate from other hemi
        if hemi == 'lh':
            otherhemi = 'rh'
        else:
            otherhemi = 'lh'
        var_x = loadtxt('%s/%s/fit_x-%s.txt'%(datadir,sub,otherhemi))
        var_y = loadtxt('%s/%s/fit_y-%s.txt'%(datadir,sub,otherhemi))
    
    # define weights per connection from variogram
    w_arr = getWeightMatrix(adj,var_x,var_y)

    # get fc value matrix
    data = loadtxt('%s/fc-labels-%s.txt'%(datadir_func,hemi))
    print(shape(data))

    # get residuals after regressing out these weights related to spatial autocorrelation
    res_arr = getResiduals(data,w_arr)

    savetxt('%s/fc-labels-%s-sareg.txt'%(datadir_func,hemi),res_arr)

    # clean up the big files
    dist_file='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/cortex/adj3_uint16-cortex-32k-%s.npz'%(sub,hemi)
    try:
        os.remove(dist_file)
    except:
        pass

    fc_file = '/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/cortex/fc_32k232k_0mm-bp/relmatch/sub-%s/ses-%s/hemi-%s/fc-%s_cortex_hybrid-%s.npz'%(sub,sub,sub[-1],hemi,hemi,hemi)
    try:
        os.remove(fc_file)
    except:
        pass

