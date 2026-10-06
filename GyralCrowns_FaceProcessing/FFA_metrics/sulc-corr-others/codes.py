import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from matplotlib import pyplot as plt


def meas_corr_meas(meas1, meas2):
    figsize = (6, 5)
    hemis = ('lh', 'rh')
    rois = ('mFus', 'pFus')
    fname1 = f'FFA_{meas1}.csv'
    fname2 = f'FFA_{meas2}.csv'
    out_name = f'{meas1}-corr-{meas2}'
    out_file = f'{out_name}.jpg'
    meas2name = {
        'va': 'surface area (mm^2)',
        'sulc': 'sulcus depth (cm)',
        'sulcBtm': 'sulcus bottom (cm)',
        'gyralCrown': 'gyral crown (cm)',
        'activ': 'face selectivity (WM)',
        'activ-emo': 'face selectivity (Emo)'}

    n_hemi = len(hemis)
    n_roi = len(rois)
    df1 = pd.read_csv(fname1)
    df2 = pd.read_csv(fname2)

    fig, axes = plt.subplots(n_roi, n_hemi, figsize=figsize)
    for roi_idx, roi in enumerate(rois):
        for hemi_idx, hemi in enumerate(hemis):
            ax = axes[roi_idx, hemi_idx]
            col = f'{hemi}_{roi}'
            x = df1[col].values
            y = df2[col].values
            nan_vec = np.logical_or(np.isnan(x), np.isnan(y))
            non_nan_vec = ~nan_vec
            x = x[non_nan_vec]
            y = y[non_nan_vec]

            print(f'{out_name} in {col}:\n',
                  pearsonr(x, y, alternative='two-sided'))
            # plot scatter
            ax.scatter(x, y, c='k', s=3)

            # fit and construct polynomial
            coefs = np.polyfit(x, y, 1)
            polynomial = np.poly1d(coefs)

            # plot fitted line
            x_min, x_max = np.min(x), np.max(x)
            x_plot = np.linspace(x_min, x_max, 100)
            y_plot = polynomial(x_plot)
            ax.plot(x_plot, y_plot, color='r')

            if roi_idx == 0:
                ax.set_title(hemi)
            elif roi_idx == n_roi - 1:
                ax.set_xlabel(meas2name[meas1])
            if hemi_idx == 0:
                ax.set_ylabel(f'{roi}\n{meas2name[meas2]}')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_file)


if __name__ == '__main__':
    meas_corr_meas(meas1='sulc', meas2='va')
    meas_corr_meas(meas1='sulc', meas2='activ')
    meas_corr_meas(meas1='sulc', meas2='activ-emo')
    meas_corr_meas(meas1='sulcBtm', meas2='va')
    meas_corr_meas(meas1='sulcBtm', meas2='activ')
    meas_corr_meas(meas1='sulcBtm', meas2='activ-emo')
    meas_corr_meas(meas1='gyralCrown', meas2='va')
    meas_corr_meas(meas1='gyralCrown', meas2='activ')
    meas_corr_meas(meas1='gyralCrown', meas2='activ-emo')
