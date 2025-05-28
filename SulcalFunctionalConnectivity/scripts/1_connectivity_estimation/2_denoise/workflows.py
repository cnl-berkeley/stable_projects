import os
from nipype.pipeline import engine as pe
from niworkflows.interfaces import confounds
import nipype.interfaces.io as nio
import nipype.interfaces.utility as niu
import nipype.interfaces.fsl as fsl
from nipype.interfaces import afni, base
import utils

from nipype import config
config.enable_debug_mode()

def activation_wf(subject, session, task, hemi, run_inds, opts, name='wf_activation'):

    print(f'\tParameters:')
    print(f'\t\tDropvols: {opts.dropvols}')
    print(f'\t\tFriston24: {opts.friston24}')
    print(f'\t\taCompCor: {opts.acompcor}')
    print(f'\t\ttCompCor: {opts.tcompcor}')
    print(f'\t\tCSF/WM: {opts.tissue_signals}')
    print(f'\t\tGlobal Signal: {opts.globalsignal}')
    print(f'\t\tHighpass: {opts.highpass} Hz')
    print(f'\t\tDetrend: {True}')
    print(f'\t\tOutlier FD threshold: {opts.outlier_def_fd}')
    print(f'\t\tOutlier DVARS threshold: {opts.outlier_def_dvars}')
    print(f'\t\tOutlier flag N scans: {opts.outlier_flag_inds}')
    print('*',run_inds)

    if not os.path.exists(opts.output_dir): os.makedirs(opts.output_dir, exist_ok=True)

    workflow = pe.Workflow(name=name)

    datasink = pe.Node(nio.DataSink(), name=f'{name}_datasink')
    if session != '':
        datasink.inputs.base_directory = '%s/sub-%s/ses-%s/hemi-%s'%(opts.output_dir,subject,session,hemi)
    else:
        datasink.inputs.base_directory = '%s/sub-%s/hemi-%s'%(opts.output_dir,subject,hemi)
    datasink.inputs.substitutions = [('_run_', 'run-0'),
                                    ('_flameo','flameo')]

    # Set up input
    inputnode = pe.Node(interface=niu.IdentityInterface(fields=['subject', 'session', 'run', 'task', 'hemi'],
                        mandatory_inputs=True), name=f'{name}_inputnode')
    inputnode.inputs.subject = subject
    inputnode.inputs.session = session
    inputnode.iterables = [('run', run_inds)]
    inputnode.inputs.task = task
    inputnode.inputs.hemi = hemi

    datasource_func = pe.Node(interface=niu.Function(function=utils.get_func_inputs,
                        input_names=['bids_dir','derivatives_dir','subject','session','run','task','hemi'],
                        output_names=['prepped_bold','confounds_file','events_file','tr']),
                        name=f'{name}_datasource-func')
    datasource_func.inputs.bids_dir = opts.bids_dir
    datasource_func.inputs.derivatives_dir = opts.derivatives_dir

    datasource_struc = pe.Node(interface=niu.Function(function=utils.get_struc_inputs,
                        input_names=['derivatives_dir','fs_subjects_dir','subject','session','hemi'],
                        output_names=['surface_file','freesurfer_dir','hemi']),
                        name=f'{name}_datasource-struc')
    datasource_struc.inputs.bids_dir = opts.bids_dir
    datasource_struc.inputs.derivatives_dir = opts.derivatives_dir
    datasource_struc.inputs.fs_subjects_dir = opts.fs_subjects_dir
    datasource_struc.inputs.hemi = hemi
    datasource_struc.inputs.subject = subject
    datasource_struc.inputs.session = session

    wb_surfsmooth = pe.Node(interface=niu.Function(function=utils.wb_surfsmooth,
                                input_names=['gii_file','surface_file','fwhm'],
                                output_names=['smooth_gii_file']),
                                name=f'{name}_surfsmooth')
    wb_surfsmooth.inputs.fwhm = opts.smooth_surf_mm

    wb_gii2nii = pe.Node(interface=niu.Function(function=utils.wb_gii2nii,
                                input_names=['gii_file'],
                                output_names=['nii_file']),
                                name=f'{name}_gii2nii')

    wb_nii2gii = pe.Node(interface=niu.Function(function=utils.wb_nii2gii,
                                input_names=['gii_file','surface_file'],
                                output_names=['nii_file']),
                                name=f'{name}_nii2gii')

    # Truncation
    # requires nifti format
    slicer = pe.Node(interface=fsl.utils.ExtractROI(), output_type='NIFTI_GZ', name=f'{name}_temporal-roi-extractor')
    slicer.inputs.t_min = opts.dropvols
    slicer.inputs.t_size = -1

    # Despiking
    # requires nifti format
    despiker = pe.Node(interface=afni.Despike(), name=f'{name}_despike')
    despiker.inputs.outputtype = 'NIFTI_GZ'
    despiker.inputs.args = '-nomask' # process all voxels

    wb_nii2gii = pe.Node(interface=niu.Function(function=utils.wb_nii2gii,
                                input_names=['nii_file', 'surface_file'],
                                output_names=['gii_file']),
                                name=f'{name}_nii2gii')

    # Estimating motion spikes
    outlier_estimator = pe.Node(interface=confounds.SpikeRegressors(),
                        name=f'{name}_outlier-estimator')
    outlier_estimator.inputs.concatenate = False
    outlier_estimator.inputs.fd_thresh = opts.outlier_def_fd
    outlier_estimator.inputs.dvars_thresh = opts.outlier_def_dvars
    if opts.outlier_flag_inds != 0:
        outlier_estimator.inputs.lags = list(range(1,opts.outlier_flag_inds+1,1)) # flag volumes after spike
    else:
        outlier_estimator.inputs.lags = [0]
    outlier_estimator.inputs.output_format = 'spikes'
    outlier_estimator.inputs.header_prefix = 'manual_outlier'

    # Setting up confounds
    confound_preprocessor = pe.Node(interface=niu.Function(function=utils.prepare_confounds_file,
                                input_names=['confound_filepath', 'spike_filepath', 'dropvols', 'acompcor', 'tcompcor', 'tissue_signals', 'highpass', 'friston24', 'global_signal'],
                                output_names=['confounds_txt', 'confounds_csv', 'outlier_inds_csv']), 
                                name=f'{name}_confound-preprocessor')
    confound_preprocessor.inputs.friston24 = opts.friston24
    confound_preprocessor.inputs.dropvols = opts.dropvols
    confound_preprocessor.inputs.tissue_signals = opts.tissue_signals
    confound_preprocessor.inputs.highpass = opts.highpass
    confound_preprocessor.inputs.acompcor = opts.acompcor
    confound_preprocessor.inputs.tcompcor = opts.tcompcor
    confound_preprocessor.inputs.global_signal = opts.globalsignal

    # Setting up bandpassing
    confound_bandpasser = pe.Node(interface=utils.OneDBandpass(), name=f'{name}_confound-bandpasser')
    confound_bandpasser.inputs.no_detrend = True
    confound_bandpasser.inputs.highpass = opts.highpass
    confound_bandpasser.inputs.lowpass = opts.lowpass

    nifti_bandpasser = pe.Node(interface=afni.Bandpass(), name=f'{name}_nifti-bandpasser')
    nifti_bandpasser.inputs.no_detrend = True
    nifti_bandpasser.inputs.outputtype = 'NIFTI_GZ'
    nifti_bandpasser.inputs.highpass = opts.highpass
    nifti_bandpasser.inputs.lowpass = opts.lowpass

    bp_mean_image_calc = pe.Node(interface=fsl.maths.MeanImage(), name=f'{name}_bandpasser-mean-image-calc')
    bp_mean_image_calc.inputs.dimension = 'T'

    remeaner = pe.Node(niu.Function(function=utils.remean_image,
                                                input_names=['demeaned_image', 'image_with_mean'],
                                                output_names=['out_file']),
                                                name=f'{name}_remeaning-node'
                                                )

    clean_confound_appender = pe.Node(interface=niu.Function(function=utils.demean_and_add_regressors,
                input_names=['txt_path', 'csv_path', 'add_linear_term', 'add_intercept_term'],
                output_names=['confounds_txt', 'confounds_csv']), name=f'{name}_clean-confound-appender')
    if opts.no_detrend:
        clean_confound_appender.inputs.add_linear_term = False
    else:
        clean_confound_appender.inputs.add_linear_term = True


    # Set up GLM
    glm_denoising = pe.Node(niu.Function(function=utils.nilearn_regress,
                                      input_names=['in_file','design'],
                                      output_names=['out_file']),
                                      name=f'{name}-nilearn-regression')

    num_copes = lambda x: len(x)

    print('\tConnecting workflow.')
    #workflow.write_graph()
    workflow.connect([
        (inputnode, datasource_func, [('subject', 'subject'),
                                ('session', 'session'),
                                ('run', 'run'),
                                ('task', 'task'),
                                ('hemi','hemi')]),
        (datasource_func, wb_surfsmooth, [('prepped_bold','gii_file')]),
        (datasource_struc, wb_surfsmooth, [('surface_file','surface_file')]),
        (wb_surfsmooth, wb_gii2nii, [('smooth_gii_file','gii_file')]),
        (wb_gii2nii, slicer, [('nii_file','in_file')])
    ])

    workflow.connect([
        (datasource_func, outlier_estimator, [('confounds_file', 'confounds_file')]),
        (outlier_estimator, confound_preprocessor, [('confounds_file', 'spike_filepath')]),
        (outlier_estimator, datasink, [('confounds_file', 'lev1.@outliers')]),
        (datasource_func, confound_preprocessor, [('confounds_file', 'confound_filepath')]),
        (confound_preprocessor, datasink, [('confounds_csv', 'lev1.@confounds_csv'),
                                           ('outlier_inds_csv','lev1.@outlier_inds_csv')]),
    ])

    if opts.lowpass != 0 or opts.highpass != 0:
        if opts.no_despike:
            workflow.connect([
                (slicer, bp_mean_image_calc, [('roi_file', 'in_file')]),
                (slicer, nifti_bandpasser, [('roi_file', 'in_file')])
            ])
        else:
            workflow.connect([
                (slicer, despiker, [('roi_file', 'in_file')]),
                (despiker, bp_mean_image_calc, [('out_file', 'in_file')]),
                (despiker, nifti_bandpasser, [('out_file', 'in_file')])
            ])

        workflow.connect([
            (datasource_func, confound_bandpasser, [('tr', 'dt')]),
            (nifti_bandpasser, remeaner, [('out_file', 'demeaned_image')]),
            (bp_mean_image_calc, remeaner, [('out_file', 'image_with_mean')]),
            (remeaner, glm_denoising, [('out_file', 'in_file')]),
            (confound_preprocessor, confound_bandpasser, [('confounds_txt', 'in_file')]),
            (confound_preprocessor, clean_confound_appender, [('confounds_csv', 'csv_path')]),
            (confound_bandpasser, clean_confound_appender, [('out_file', 'txt_path')]),
            (clean_confound_appender, glm_denoising, [('confounds_txt', 'design')]),
            (datasource_struc, wb_nii2gii, [('surface_file','surface_file')]),
            (glm_denoising, wb_nii2gii, [('out_file','nii_file')]),
            (wb_nii2gii, datasink, [('gii_file','lev1.@residuals_gii')])
        ])

    else:
        sys.exit('Not-filtering not implemented for this pipeline.')

    return workflow

