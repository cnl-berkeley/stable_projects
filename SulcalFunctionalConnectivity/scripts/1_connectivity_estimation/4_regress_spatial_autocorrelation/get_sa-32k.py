import os,sys
from numpy import *
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


def getInputs(sub,hemi,outdir):
    dist_file='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/cortex/adj3_uint16-cortex-32k-%s.npz'%(sub,hemi)
    arr = load(dist_file,allow_pickle=True)
    adj = arr['arr_0']
    dist_vec = adj.flatten()

    print(shape(adj),shape(dist_vec))

    # get fc value matrix
    datadir='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/cortex/fc_32k232k_0mm-bp-relmatch'%sub
    fc_file = '%s/sub-%s/ses-%s/hemi-%s/fc-%s_cortex_hybrid-%s.npz'%(datadir,sub,sub[-1],hemi,hemi,hemi)
    arr = load(fc_file,allow_pickle=True)
    fc_arr = arr['arr_0']
    fc_vec = fc_arr.flatten()
    print(shape(fc_arr),shape(fc_vec))

    # exclude long-range connections (to save memory)
    i1 = argwhere(dist_vec > 302) # 3 cm

    # exclude self-connections    
    i2 = argwhere(dist_vec == 0) # 32k has ca 2 mm spacing
    i_ex = concatenate((i1,i2))
    print(shape(i_ex))

    fc_vec = delete(fc_vec,i_ex)
    dist_vec = delete(dist_vec,i_ex)

    print(shape(fc_vec),amin(fc_vec),amax(fc_vec))
    print(shape(dist_vec),amin(dist_vec),amax(dist_vec))

    return fc_vec, dist_vec


def plotRaw(outdir,hemi,dist_vec,fc_vec):
    import pylab as P
    import matplotlib.pyplot as plt
    font0 = FontProperties(size=4)
    fig = P.figure(figsize=(12,12))
    x1 = fig.add_axes([0.25, 0.25, 0.46, 0.23])

    plt.plot(dist_vec,fc_vec, '.b',alpha=0.1,markersize=1)

    def myformat(x, pos):
        return x

    x1.xaxis.set_major_locator(MultipleLocator(10))
    x1.xaxis.set_major_formatter(myformat)
    x1.xaxis.set_minor_locator(MultipleLocator(1))

    x1.tick_params(axis='both', which='major', length=5)

    fig.savefig('%s/raw-32k-%s.png'%(outdir,hemi),bbox_inches='tight')
    P.close(fig)


def spatial_autocorrelation(dist_flat,cm_flat, discretization=5):
    # https://github.com/murraylab/spatiotemporal/blob/main/spatiotemporal/stats.py
    import pandas
    import numpy as np
    import scipy.spatial
    import scipy.optimize

    df = pandas.DataFrame(np.asarray([dist_flat, cm_flat]).T, columns=["dist", "corr"])

    df['dist_bin'] = np.round(df['dist']/discretization)*discretization
    bin_means = df.groupby('dist_bin').mean().values
    print('bin_means',bin_means)
    df_binned = df.groupby('dist_bin').mean().reset_index().sort_values('dist_bin')
    binned_dist_flat = df_binned['dist_bin']
    print(shape(binned_dist_flat),binned_dist_flat)
    binned_cm_flat = df_binned['corr']
    print(shape(binned_cm_flat))
    spatialfunc = lambda v : exp((-binned_dist_flat+v[2])/v[0])*(1-v[1])+v[1]
    with np.errstate(all='warn'):
        res = scipy.optimize.minimize(lambda v : sum((binned_cm_flat-spatialfunc(v))**2), [30, .1, 10], bounds=[(.1, 100), (-1, 1), (-20, 20)])

    print(res)

    return (res.x[0],res.x[1],res.x[2],bin_means,binned_dist_flat,binned_cm_flat)


if __name__ == '__main__':

    sub = sys.argv[1]
    hemi = sys.argv[2]

    outdir='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/cortex'%sub
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    fc_vec,dist_vec = getInputs(sub,hemi,outdir)

    plotRaw(outdir,hemi,dist_vec,fc_vec)

    v = spatial_autocorrelation(dist_vec,fc_vec)
    print(v[0],v[1],v[2])

    xdata = dist_vec.flatten()
    ydata = fc_vec.flatten()

    xi = v[3]
    yi = exp(-1.0*(xi-v[2])/v[0])*(1-v[1])+v[1]

    xi2 = arange(0,300,2)
    yi2 = exp(-1.0*(xi2-v[2])/v[0])*(1-v[1])+v[1]

    print('xi',shape(xi),xi)
    print('yi',shape(yi),yi)

    savetxt('%s/fit_x-%s.txt'%(outdir,hemi),xi2)
    savetxt('%s/fit_y-%s.txt'%(outdir,hemi),yi2)

    import pylab as P
    import matplotlib.pyplot as plt
    fig = P.figure(figsize=(12,12))
    x1 = fig.add_axes([0.25, 0.25, 0.46, 0.23])

    plt.plot(xdata,ydata,'.g',alpha=0.2,markersize=0.5)
    plt.plot(xi2, yi2, '-k',lw=1,alpha=1)

    def myformat(x, pos):
        y = 0.1*x
        return y.astype(int)

    x1.xaxis.set_major_locator(MultipleLocator(100))
    x1.xaxis.set_major_formatter(myformat)
    x1.xaxis.set_minor_locator(MultipleLocator(10))

    x1.tick_params(axis='both', which='major', length=5)

    x1.spines['left'].set_linewidth(0.5)
    x1.spines['right'].set_linewidth(0)
    x1.spines['top'].set_linewidth(0)
    x1.spines['bottom'].set_linewidth(0.5)

    fig.savefig('%s/fit-32k-%s2.png'%(outdir,hemi),bbox_inches='tight')
