import nipype.interfaces.workbench.base as wb
import nibabel as nib
from numpy import *
import os,sys

os.environ['SUBJECTS_DIR'] = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7'

task = 'relmatch'
subject_file='/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/scripts/fmri_cohort/relmatch_parlabels3/namelist_full.txt'
subs = loadtxt(subject_file,dtype='str')

def incCiter(citer):
    if citer[0] < 255:
        citer[0] = citer[0]+1
    elif citer[1] < 255:
        citer[1] = citer[1]+1
    elif citer[2] < 255:
        citer[2] = citer[2]+1
    else:
        sys.exit('No color available.')
    return citer

# population probabilistic maps (MIDB Precision Brain Atlas)
atlasname = 'combined_clusters_thresh0.61'
dlabelpath = '%s/%s.dlabel.nii'%(os.getcwd(),atlasname)

# separate to lh and rh cortex
annot_lh = '%s/lh.%s.label.gii'%(os.getcwd(),atlasname)
annot_rh = '%s/rh.%s.label.gii'%(os.getcwd(),atlasname)
str = 'wb_command -cifti-separate %s COLUMN -label CORTEX_LEFT %s -label CORTEX_RIGHT %s'%(dlabelpath,annot_lh,annot_rh)
cmd=wb.WBCommand(str)
cmd.run()

# parse names (e.g. for network identity)
labelnames = '%s/%s.names.csv'%(os.getcwd(),atlasname)
str = 'wb_command -cifti-label-export-table %s 1 %s'%(dlabelpath,labelnames)
cmd=wb.WBCommand(str)
cmd.run()

# create freesurfer format annotation files
inputf = open(labelnames)
lines = inputf.readlines()
inputf.close()

for hemi in ['lh','rh']:
    im_file=nib.load(annot_lh)
    ts=squeeze(asarray([x.data for x in im_file.darrays]))
    n_vinds = len(ts)

    label_array = array([0]*n_vinds)
    ctab = []
    citer = [2,0,0,255]
    networks = []
    for i in arange(0,len(lines),2):
        tmp = lines[i+1].split()
        labelno = int(tmp[0])
        color = citer
        citer = incCiter(citer)

        if i == 0: # with first label, add also 0 unknown
            ctab.append([1,0,0,0])
            names = ['0-%s.unknown'%hemi]

        label_array[ts == labelno] = labelno
        ctab.append(citer.copy())

        tmp = '%d-%s.%s'%(labelno,hemi,lines[i].strip('\n'))
        names.append(tmp)

    nw_ctab = []
    nw_citer = [2,0,0,255]
    uni_labels = unique(label_array)
    label_array_nw = array([0]*n_vinds)
    nw_names = []
    for i in arange(len(uni_labels)):
        nw_name = names[i].split('.')[-1].split('_')[0]
        print(i,names[i],nw_name)
        try:
            nw_ind = nw_names.index(nw_name)
        except:
            nw_names.append(nw_name)
            nw_ind = len(nw_names)-1
            nw_ctab.append(nw_citer.copy())
            nw_citer = incCiter(nw_citer)
        
        label_array_nw[label_array == uni_labels[i]] = nw_ind 

    print(len(unique(label_array_nw)),len(unique(nw_names)))
    filepath = '%s/%s.%s-nw-cifti32k.annot'%(os.getcwd(),atlasname,hemi)
    nib.freesurfer.io.write_annot(filepath, label_array_nw, array(nw_ctab), nw_names, fill_ctab=True)

    # since the atlas is defined on the standard 32k surface, convert to full surface (native)

    import nipype.interfaces.freesurfer as fs
    lowres_sphere = '%s/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.%s.sphere.32k_fs_LR.surf.gii'%(os.getcwd(),hemi[0].upper())
    sphere_bin = '%s/%s.sphere32k'%(os.getcwd(),hemi)
    mris = fs.MRIsConvert()
    mris.inputs.in_file = lowres_sphere
    mris.inputs.out_file = sphere_bin
    mris.run()

    surf_data = nib.freesurfer.io.read_geometry(sphere_bin)
    surf_vcoords = surf_data[0]
    surf_vinds = arange(0,max(surf_data[1].flatten())+1,1)

    # loop subjects
    for sub in subs:
        full_surf_file = '%s/sub-%s/surf/%s.sphere'%(os.environ['SUBJECTS_DIR'],sub,hemi)
        full_surf_data = nib.freesurfer.io.read_geometry(full_surf_file)
        full_surf_vcoords = full_surf_data[0]
        full_surf_vinds = arange(0,max(full_surf_data[1].flatten())+1,1)

        annot_32k = '%s/%s.%s-nw-cifti32k.annot'%(os.getcwd(),atlasname,hemi)
        labels1,ctab1,names1 = nib.freesurfer.io.read_annot(annot_32k)
        print(len(unique(labels1)),len(unique(names1)))

        label_array = array([0]*len(full_surf_vcoords))

        for i in arange(len(full_surf_vcoords)):
            c1 = full_surf_vcoords[i]
            dist_arr = sqrt((c1[0]-surf_vcoords[:,0])**2 + (c1[1]-surf_vcoords[:,1])**2 + (c1[2]-surf_vcoords[:,2])**2)
            label_array[i] = labels1[argmin(dist_arr)]

            # annot colortable and names remain the same as in 32k resolution
        savedir = '%s/annot-%s'%(os.getcwd(),task)
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        filepath = '%s/sub-%s_ses-%s_%s-%s-nw-full.annot'%(savedir,sub,sub[-1],atlasname,hemi)
        annot = nib.freesurfer.io.write_annot(filepath, label_array, array(ctab1), names1, fill_ctab=True)

