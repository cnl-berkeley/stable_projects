import sys
import argparse
import os
import datetime
import subprocess
import pandas as pd
from numpy import arange

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('input_directory',default='/home/weiner/Nora_PFCSulci/BIDS/derivatives/fmriprep')
    parser.add_argument('--output_directory',default='/home/weiner/shakki/fmriprep/derivatives')
    parser.add_argument('--logdir')
    parser.add_argument('--jobname')
    parser.add_argument('--idfile')
    parser.add_argument('--id')
    parser.add_argument('--qsub-template')

    parser.add_argument('--fs_lic_path',default='/usr/local/freesurfer_x86_64-7.1.0/license.txt')
    parser.add_argument('--fs_subjects_dir',default='/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/fs_v7_16-18')

    parser.add_argument('--mem-gb', type=int, default=2)
    parser.add_argument('--num-cores', type=int, default=2)
    parser.add_argument('--num-hours', type=int, default=10)

    args = parser.parse_args()
    if os.path.exists(args.input_directory) is False:
        print('Input directory {} does not exist'.format(args.input_directory))
        sys.exit(1)

    if args.output_directory is None:
        output_directory = os.path.join(os.path.dirname(os.path.dirname(args.input_directory)), 'fmriprep')
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)
        args.output_directory = output_directory

    if args.logdir is None:
        logdir = os.path.join(os.path.dirname(os.path.dirname(args.input_directory)), 'logs')
        if os.path.exists(logdir) is False:
            os.mkdir(logdir)
        args.logdir = logdir

    if args.logdir and os.path.exists(args.logdir) is False:
        os.mkdir(args.logdir)
    
    if args.jobname is None:
        args.jobname = 'fmriprep-job'

    return args

def parse_idfile(idfile):
    data=pd.read_csv(idfile,header=None,sep='t',names=['subject','session'])
    data['subject'] = data['subject'].apply(lambda x: x[1::]) # remove n from front of identifier
    folder_run_list=[]
    for i in arange(len(data['subject'])):
        sub='sub-%s'%data['subject'][i]
        session='ses-%s'%data['session'][i]
        folder_run_list.append([sub,session])
    return folder_run_list

def write_run_list(run_list):
    now = datetime.datetime.now()
    filename = '{}{}{}-{}{}-{}_fmriprep_run_files-list'.format(now.year, str(now.month).zfill(2), str(now.day).zfill(2), now.hour,  now.minute, now.second, now.microsecond)
    if os.path.exists(filename): os.remove(filename)
    with open(filename, 'w') as f:
        for itm in run_list:
            f.write(' '.join(itm) + '\n')
    return os.path.abspath(filename)

def main():
    args = parse_args()

    folder_run_list = parse_idfile(args.idfile)
    run_list_abspath = write_run_list(folder_run_list)

    runvars = {'num_jobs':    len(folder_run_list),
               'num_cores':   args.num_cores,
               'num_threads': args.num_cores * 2,
               'num_hours':   '{}:00:00'.format(args.num_hours),
               'mem_gb':      args.mem_gb,
               'mem_mb':      args.mem_gb * 1000 * args.num_cores,
               'logdir':      args.logdir,
               'jobname':     args.jobname,
               'runlist':     run_list_abspath,
               'root':        args.input_directory,
               'outputpath':  args.output_directory,
               'fs_lic_dir':  args.fs_lic_dir,
               'fs_subjects_dir': args.fs_subjects_dir
               }
    
    print(runvars)
    cmd = 'qsub -t 1-{num_jobs} -S /bin/bash -j y -o {logdir} -e {logdir} -pe smp {num_cores} -l mem_free={mem_gb} -l h_rt={num_hours} -N {jobname} -v MEM_IN_MB={mem_mb} -v LOGDIR={logdir} -v NUM_THREADS={num_threads} -v OUTPUTDATAPATH={outputpath} -v BIDS_DATA_DIRECTORY={root} -v FS_LIC_PATH={fs_lic_path} -v FS_SUBJECTS_DIR={fs_subjects_dir} -v RUNLIST={runlist} {script}'.format(**runvars)
    subprocess.call(cmd.split())

if __name__ == "__main__":
    main()

