import nibabel as nib
import numpy as np
import pandas as pd
import traits
from nipype.interfaces import afni, base

def get_func_inputs(bids_dir,derivatives_dir,subject,session,run,task,hemi):
    import sys
    from glob import glob
    import json

    if session != '':
        sub_derivatives_dir = '%s/sub-%s/ses-%s'%(derivatives_dir,subject,session)
    else:
        sub_derivatives_dir = '%s/sub-%s'%(derivatives_dir,subject)

    hemiletter = hemi[0].upper()
    prepped_bold = glob('%s/func/*task-%s_run-%d_hemi-%s_space-fsnative_bold.func.gii'%(sub_derivatives_dir,task,run,hemiletter))
    prepped_json = glob('%s/func/*task-%s_run-%d_hemi-%s_space-fsnative_bold.json'%(sub_derivatives_dir,task,run,hemiletter))
    if not prepped_bold or not prepped_json:
        sys.exit('No preprocessed files found.')
    with open(prepped_json[0],'r') as json_file:
        json_content = json.load(json_file)
        tr = json_content['RepetitionTime']
    print('tr',tr)

    confounds_file = glob('%s/func/*task-%s_run-%d_desc-confounds_timeseries.tsv'%(sub_derivatives_dir,task,run))
    if not confounds_file:
        sys.exit('No confounds file found.')

    events_file = glob('%s/*task-%s_run-%d_events.tsv'%(bids_dir,task,run))
    if not events_file:
        sys.exit('No events file found.')

    return prepped_bold[0],confounds_file[0],events_file[0],tr


def get_struc_inputs(derivatives_dir,fs_subjects_dir,subject,session,hemi):
    import sys
    from glob import glob

    if session != '':
        sub_derivatives_dir = '%s/sub-%s/ses-%s'%(derivatives_dir,subject,session)
        surface_file = glob('%s/anat/sub-%s_ses-%s_hemi-%s_midthickness.surf.gii'%(sub_derivatives_dir,subject,session,hemi[0].upper()))
    else:
        sub_derivatives_dir = '%s/sub-%s'%(derivatives_dir,subject)
        surface_file = glob('%s/anat/sub-%s_hemi-%s_midthickness.surf.gii'%(sub_derivatives_dir,subject,hemi[0].upper()))
    if not surface_file:
        sys.exit('No surface file found.')

    freesurfer_dir = '%s/sub-%s'%(fs_subjects_dir,subject)

    return surface_file[0],freesurfer_dir,hemi


def prepare_confounds_file(confound_filepath, spike_filepath, dropvols, acompcor, tcompcor, tissue_signals, highpass, friston24, global_signal):
    import os
    import pandas as pd
    import numpy as np

    confounds_df = pd.read_csv(confound_filepath, sep='\t')
    confounds_df = confounds_df.fillna(0)
    confounds_to_use = []

    friston24_cols = [c for c in confounds_df.columns if 'rot' in c.lower() or 'trans' in c.lower()]
    if friston24:
        print('Using Friston 24 as motion regressors.')
        confounds_to_use += friston24_cols
    if (acompcor or tcompcor) and highpass == 0:
        print('Including cosines (using CompCors without other high-pass filtering).')
        confounds_to_use += [c for c in confounds_df.columns if 'cosine' in c.lower()]
    if acompcor:
        print('Using aCompCor.')
        all_acompcor_cols = [c for c in confounds_df.columns if 'a_comp' in  c.lower()]
        all_acompcor_cols = sorted(all_acompcor_cols, key=lambda colname: int(colname.split('_')[-1]))
        confounds_to_use += all_acompcor_cols[:5] # first 5 compcor
    if tcompcor:
        print('Using tCompCor.')
        confounds_to_use += [c for c in confounds_df.columns if 't_comp' in c.lower()]
    if tissue_signals:
        print('Using mean tissue signals from eroded segmentation.')
        wm_csf_signals = [c for c in confounds_df.columns if 'csf' in c.lower() or 'white_matter' in c.lower()]
        wm_csf_signals = [c for c in wm_csf_signals if 'square' not in c and 'deriv' not in c and 'power' not in c]
        confounds_to_use += wm_csf_signals
    if global_signal:
        print('Including global signal.')
        confounds_to_use.append('global_signal')
    confounds_filtered = confounds_df[confounds_to_use][dropvols:]
    confounds_filtered = confounds_filtered.round(decimals=3)

    # truncate and format info on motion spikes
    if os.path.exists(spike_filepath):
        try:
            spike_df = pd.read_csv(spike_filepath, sep='\t')
            spike_cols = [c for c in spike_df.columns if 'manual_outlier' in c.lower()]
            spike_df = spike_df[spike_cols]
            outliermask = spike_df.sum(1)
            outliermask = np.abs(outliermask[dropvols:]-1)
            outlier_inds = [i for i, d in enumerate(outliermask) if d == 0]
        except: # file not found or no spikes
            outlier_inds = pd.DataFrame()

    else:
        outlier_inds = pd.DataFrame()

    outlier_inds_csv = '%s/event-inds_outliers.csv'%os.getcwd() # zero-based indices
    np.savetxt(outlier_inds_csv, outlier_inds, delimiter=', ', fmt='%d')

    filepath_split = os.path.basename(confound_filepath).split('.')
    filepath_stem = ''.join(filepath_split[:-1]).split('desc')
    filepath_stem = ''.join(filepath_stem[:1]) + 'desc'

    out_txtfile_path = os.path.join(os.getcwd(), f'{filepath_stem}_confounds.txt')
    out_csvfile_path = os.path.join(os.getcwd(), f'{filepath_stem}_confounds.csv')
    
    np.savetxt(out_txtfile_path, confounds_filtered.values)
    confounds_filtered.to_csv(out_csvfile_path, index=False)

    return out_txtfile_path,out_csvfile_path,outlier_inds_csv


def demean_and_add_regressors(txt_path, csv_path, add_linear_term=True, add_intercept_term=False):
    from numpy import genfromtxt, savetxt, c_, arange, ones
    import os
    from pandas import read_csv

    txt_array = genfromtxt(txt_path)
    csv_as_df = read_csv(csv_path)
    demeaned_txt_array = txt_array + csv_as_df.values.mean(axis=0)

    csv_as_df = read_csv(csv_path)
    csv_as_df[:] = demeaned_txt_array

    demeaned_csv_as_df = csv_as_df

    num_rows = len(txt_array)
    if add_linear_term:
        linear_term = arange(1, num_rows+1)
        demeaned_txt_array = c_[demeaned_txt_array, linear_term]
        demeaned_csv_as_df['linear_trend'] = linear_term
        print(f'Appending linear term to design. Shape of design matrix was: {txt_array.shape}, now {demeaned_txt_array.shape}.')

    if add_intercept_term:
        intercept_term = ones(num_rows)
        demeaned_txt_array = c_[demeaned_txt_array, intercept_term]
        demeaned_csv_as_df['intercept_term'] = intercept_term
        print(f'Appending intercept term to design. Shape of design matrix was: {txt_array.shape}, now {demeaned_txt_array.shape}.')

    filepath_split = os.path.basename(csv_path).split('.')
    filepath_stem = ''.join(filepath_split[:-1]).split('desc')
    filepath_stem = ''.join(filepath_stem[:1]) + 'desc'

    out_txtfile_path = os.path.join(os.getcwd(), f'{filepath_stem}-cleanconfounds_bold.txt')
    out_csvfile_path = os.path.join(os.getcwd(), f'{filepath_stem}-cleanconfounds_bold.csv')

    savetxt(out_txtfile_path, demeaned_txt_array)
    demeaned_csv_as_df.to_csv(out_csvfile_path, index=False)
    return out_txtfile_path, out_csvfile_path


def remean_image(demeaned_image, image_with_mean):
    from nibabel import load, Nifti1Image
    import os
    from numpy import expand_dims

    demeaned_data_nifti = load(demeaned_image)
    demeaned_data = demeaned_data_nifti.get_fdata()

    pre_unmean_nifti = load(image_with_mean)
    mean_image = pre_unmean_nifti.get_fdata()
    mean_image = expand_dims(mean_image,axis=2)

    remeaned_data = demeaned_data + mean_image[..., None]
    demeaned_image_basename = os.path.basename(demeaned_image)
    demeaned_image_base_stem = demeaned_image_basename.replace('.nii.gz', '').replace('.nii', '')

    output_basename = f'{demeaned_image_base_stem}_remeaned.nii.gz'
    output_fpath = os.path.join(os.getcwd(), output_basename)

    output_nii_object = Nifti1Image(remeaned_data, pre_unmean_nifti.affine, pre_unmean_nifti.header)
    output_nii_object.to_filename(output_fpath)

    return output_fpath


class BP1D_InputSpec(afni.base.AFNICommandInputSpec):
    in_file = base.File(
        desc='input file to bandpass',
        argstr='%s',
        position=-2,
        mandatory=True,
        exists=True,
        copyfile=True
    )

    lowpass = base.traits.Float(
        desc='lowpass',
        argstr='%f',
        position= -3,
        mandatory=True
    )

    highpass = base.traits.Float(
        desc='highpass',
        argstr='%f',
        position= -4,
        mandatory=True
    )

    dt = base.traits.Float(
        desc='TR',
        argstr='-dt %f',
        position=-5,
        mandatory=True
    )

    no_detrend = base.traits.Bool(
        argstr = '-nodetrend',
        mandatory=True,
        position=-6,
        desc='Skip the quadratic detrending that occurs before the FFT-based bandpassing'
    )

    out_file = base.File(
        desc = 'Output file from 1dBandpass',
        argstr = '> %s', 
        position= -1,
        name_source='in_file', 
        name_template='%s_bp',
        keep_extension = 'False',
    )

class BP1D_OutputSpec(base.TraitedSpec):
    out_file = base.File(exists=True, desc='output_file')

class OneDBandpass(afni.base.AFNICommand):
    _cmd = '1dBandpass'
    input_spec = BP1D_InputSpec
    output_spec = BP1D_OutputSpec


def wb_surfsmooth(gii_file,surface_file,fwhm):
    # smoothes gifti surfaces using connectome workbench
    import os
    import nipype.interfaces.workbench.base as wb

    if fwhm == 0:
        smooth_gii_file = gii_file
    else:
        smooth_gii_file = '%s/%s-smooth%dmm.func.gii'%(os.getcwd(),os.path.basename(gii_file).split('.')[0],fwhm)
        str = 'wb_command -metric-smoothing %s %s %d %s -fwhm '%(surface_file,gii_file,fwhm,smooth_gii_file)
        cmd=wb.WBCommand(str)
        cmd.run()

    return smooth_gii_file


def nilearn_regress(in_file,design):
    import nilearn.image
    
    out_data = nilearn.image.clean_img(in_file,confounds=design)

    import nibabel as nib
    import os
    out_file = '%s/res4d.nii.gz'%os.getcwd()
    nib.save(out_data,out_file)

    return out_file


def wb_gii2nii(gii_file):
    # converts metric gifti files to volume-encoded nifti format
    import os
    import nipype.interfaces.workbench.base as wb
    nii_file = '%s/%s.nii.gz'%(os.getcwd(),os.path.basename(gii_file).split('.')[0])
    str = 'wb_command -metric-convert -to-nifti %s %s'%(gii_file,nii_file)
    cmd=wb.WBCommand(str)
    cmd.run()

    return nii_file


def wb_nii2gii(nii_file,surface_file):
    import os
    import nipype.interfaces.workbench.base as wb
    gii_file = '%s/%s.func.gii'%(os.getcwd(),os.path.basename(nii_file).split('.')[0])
    str = 'wb_command -metric-convert -from-nifti %s %s %s'%(nii_file,surface_file,gii_file)
    cmd=wb.WBCommand(str)
    cmd.run()

    return gii_file

