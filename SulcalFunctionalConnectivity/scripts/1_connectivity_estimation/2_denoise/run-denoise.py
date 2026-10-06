"""

Postprocessing/denoising of data preprocessed with fMRIprep. Participant-level.

Dependencies:
# freesurfer
# fsl
# connectome workbench

# python:
# pip install niflow-nipype1-workflows
# pip install niworkflows
# bibabel

"""

#!/usr/bin/env python3
import sys
import logging
from pathlib import Path
from workflows import activation_wf
from datetime import datetime
from numpy import array

__version__ = '1.0.0'

def get_parser():
    """Define the command line interface"""
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter

    parser = ArgumentParser(description='Post-fMRIprep analysis workflow for activation analysis',
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        'derivatives_dir', action='store', type=Path,
        help='the root folder of a derivatives set generated with fMRIPrep '
             '(sub-XXXXX folders should be found at the top level in this folder).')
    parser.add_argument('output_dir', action='store', type=Path,
                        help='the output path for the outcomes of preprocessing and visual reports')

    parser.add_argument('--version', action='version', version=__version__)

    g_bids = parser.add_argument_group('Options for filtering BIDS queries')
    g_bids.add_argument('--participant-label', action='store', type=str,
                        nargs='*', help='process only particular subjects')
    g_bids.add_argument('--task', action='store', type=str, nargs='*',
                        help='select a specific task to be processed')
    g_bids.add_argument('--runstr', action='store', type=str, nargs='*',
                        help='select specific runs to be processed (list separated by commas)')
    g_bids.add_argument('--bids-dir', action='store', type=Path,
                        help='point to the BIDS root of the dataset from which the derivatives '
                             'were calculated (in case the derivatives folder is not the default '
                             '(i.e. ``BIDS_root/derivatives``).')
    g_bids.add_argument('--fs_subjects_dir', action='store', type=Path, help='FreeSurfer subject directory.')

    d_opts = parser.add_argument_group('Options for denoising.')
    d_opts.add_argument('--dropvols', type=int, default=3, help='Number of volumes to drop from the scan. Default is 3 volumes. Should match fMRIprep.')
    d_opts.add_argument('--no_detrend', action='store_true', default=False, help='Use if you do not want to do detrending.')
    d_opts.add_argument('--no_despike', action='store_true', default=False, help='Use if you do not want to do despiking using AFNI 3Despike.')
    d_opts.add_argument('--smooth_surf_mm', type=int, default=0, help='Size in mm of kernel used to smooth output data. Set to 0 if you do not want to smooth (FILM still downsamples to 2mm resolution).')
    d_opts.add_argument('--highpass', type=float, default=0.008, help='High pass filter for bandpassing. Set 0 if you do not want to high-pass filter.')
    d_opts.add_argument('--lowpass', type=float, default=0.1, help='Low pass filter for bandpassing. Set 0 if you do not want to low-pass filter.')    
    d_opts.add_argument('--acompcor', action='store_true', default=False, help='Use aCompCor. If no other high-pass filtering is specified, will also include fMRIprep cosines.')
    d_opts.add_argument('--tcompcor', action='store_true', default=False, help='Use tCompCor.')
    d_opts.add_argument('--tissue_signals', action='store_true', default=False, help='Use mean WM and CSF signals.')
    d_opts.add_argument('--friston24', action='store_true', default=False, help='Use Friston 24 head motion parameters (x, y, z displacement in mm and associated squares, derivatives of squares, and derivatives).')
    d_opts.add_argument('--globalsignal', action='store_true', default=False, help='Use to use global signal regression.')
    d_opts.add_argument('--outlier_def_fd', type=float, default=0.5, help='Minimum framewise displacement threshold to define motion outliers.')
    d_opts.add_argument('--outlier_def_dvars', type=float, default=1.5, help='Minimum DVARS threshold to define motion outliers.')
    d_opts.add_argument('--outlier_flag_inds', type=int, default=0, help='Minimum framewise displacement threshold to define motion outliers.')

    g_perfm = parser.add_argument_group('Options to handle performance')
    g_perfm.add_argument("-v", "--verbose", dest="verbose_count", action="count", default=0,
                         help="increases log verbosity for each occurence, debug level is -vvv")
    g_perfm.add_argument('--ncpus', '--nprocs', action='store', type=int,
                         help='maximum number of threads across all processes')
    g_perfm.add_argument('--nthreads', '--omp-nthreads', action='store', type=int,
                         help='maximum number of threads per-process')

    g_other = parser.add_argument_group('Other options')
    g_other.add_argument('-w', '--work-dir', action='store', type=Path,
                         help='path where intermediate results should be stored')

    return parser


def main():
#    from os import cpu_count
    from multiprocessing import set_start_method
    set_start_method('forkserver')

    opts = get_parser().parse_args()

    # Resource management options
    plugin_settings = {
        'plugin': 'MultiProc',
        'plugin_args': {
            'n_procs': 3, #opts.ncpus, set to the number of runs (Julie asked)!
            ''
            'raise_insufficient': False,
            'maxtasksperchild': 1,
        }
    }

    subject = opts.participant_label
    subject = subject[0]
    session = subject.split('t')[-1] # specify '' if there BIDS session identifier is not used
    task = opts.task[0]

    tmp = opts.runstr
    tmp = tmp[0].split(',')
    run_inds = list(array(tmp).astype(int)+1) # from 0-based to 1-based indexing
    # task conditions and contrasts are defined in utils.py in function prepare_timing_info

    for hemi in ['lh','rh']:
        print(subject, session, task, hemi, run_inds, opts)
        workflow = activation_wf(subject, session, task, hemi, run_inds, opts)
        workflow.base_dir = '%s/%s-%s-%f'%(opts.work_dir,subject,hemi,datetime.timestamp(datetime.now()))
        workflow.run(**plugin_settings)

    return 0


if __name__ == '__main__':
    sys.exit(main())

