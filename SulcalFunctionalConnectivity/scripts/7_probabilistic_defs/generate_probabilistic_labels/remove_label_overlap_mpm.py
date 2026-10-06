import os,sys
import numpy as np
import seaborn as sns
import nibabel as nib
import glob
import shutil

def getLabelFile(label,sub,hemi):
    labelfile_found = ''
    labelname_found = []
    if type(label) is list: # pick the first label in list that exits
        label_found = False
        for lab in label:
            labelfile = '%s/%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,lab)
            if os.path.exists(labelfile) and not label_found:
                label_found = True
                labelfile_found = labelfile
        labelname_found = label[0]
    else: # look for specific label
        labelfile = '%s/%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,label)
        if os.path.exists(labelfile):
            labelfile_found = labelfile
        labelname_found = label
 
    return labelfile_found,labelname_found


def getLUT(labels,label_inds_present):
    colors = sns.color_palette("hls", len(labels))
    colors = np.multiply(colors,255).astype(int)
    colors = colors[label_inds_present]
    labels = list(np.array(labels)[label_inds_present])

    lut = '%s/label_LUT.txt'%os.getcwd()
    of = open(lut,'w')
    of.write('#No. Label Name:                            R   G   B   A\n')
    for l in np.arange(len(labels)):
        print(labels[l],colors[l])
        labelname = labels[l].replace('_prob33','')
        of.write('%d   %s                                 %d   %d   %d   0\n'%(l+1,labelname,colors[l][0],colors[l][1],colors[l][2]))
    of.close()

    return lut


def label2Annot(sub,hemi,labelnames,label_files,annot_name):

    label_inds_present = []
    for i in np.arange(len(label_files)):
        try:
            if len(nib.freesurfer.read_label(label_files[i])) > 0:
                label_inds_present.append(i)
        except:
            pass

    annot_file = '%s/%s/label/%s.%s.annot'%(os.environ['SUBJECTS_DIR'],sub,hemi,annot_name)
    if os.path.exists(annot_file): # remove old
       os.remove(annot_file)

    from nipype.interfaces.freesurfer import Label2Annot
    l2a = Label2Annot()
    l2a.inputs.keep_max = True # !!!
    l2a.inputs.hemisphere = hemi
    l2a.inputs.subject_id = sub
    l2a.inputs.in_labels = list(np.array(label_files)[label_inds_present])
    l2a.inputs.orig = '%s/%s/surf/%s.pial'%(os.environ['SUBJECTS_DIR'],sub,hemi)
    l2a.inputs.out_annot = annot_name
    l2a.inputs.color_table = getLUT(labelnames,label_inds_present)
    os.system('%s --sd /home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7'%l2a.cmdline)

    return annot_file


def annot2label(sub,hemi,annot_file):
    outdir = '%s/%s/label/sulcal_prob_mpm'%(os.environ['SUBJECTS_DIR'],sub)
    cmd = 'mri_annotation2label --s %s --hemi %s --annotation %s --outdir %s'%(sub,hemi,annot_file,outdir)
    os.system(cmd)


def checkForOverlap(label_files):
    print('Checking overlap for %d labels'%len(label_files))
    n_labels = len(label_files)

    all_vertices = []
    for l1 in np.arange(n_labels):
        all_vertices.append(nib.freesurfer.read_label(label_files[l1]))    

    overlaps = np.zeros((n_labels,n_labels))*np.nan
    for l1 in np.arange(n_labels):
        for l2 in np.arange(n_labels):
            overlaps[l1,l2] = len(np.intersect1d(all_vertices[l1],all_vertices[l2]))

    print(overlaps)
    print('Vertices overlap between %d labels'%(len(np.argwhere(overlaps > 0))-n_labels))


if __name__ == '__main__':

    os.environ['SUBJECTS_DIR'] = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7'

    sub = sys.argv[1] # subject identifier
    prob_tag = sys.argv[2] # probabilistic or manually defined

    if prob_tag == 'prob':
        lpfc_labels = ['ifs','painfs_any','pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p','prts','lfms','aalf']
        lpar_labels = ['slos1','sB','pips','mTOS','iTOS','IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
        labels = lpfc_labels + lpar_labels
        for i in np.arange(len(labels)):
            labels[i] = '%s_prob33'%labels[i]
        annot_name = 'LPFC+LPC_probtest_mpm'

        outdir = '%s/%s/label/sulcal_prob_mpm'%(os.environ['SUBJECTS_DIR'],sub)
        if os.path.exists(outdir):
            shutil.rmtree(outdir)

    else:
        lpfc_labels = ['ifs',['painfs_any','painfs_combined'],'pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
        lpar_labels = [['slos1','slocs-v','SLOS'],'sB','pips','mTOS',['iTOS','ITOS','lTOS'],'IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
        mpar_labels = []
        labels = lpfc_labels + lpar_labels + mpar_labels
        annot_name = 'LPFC+LPC_mpm'

    if sub != 'fsaverage_prob':
        sub = 'sub-%s'%sub

    for hemi in ['lh','rh']:
        labelfiles = []
        labelnames = []
        for label in labels:
            labelfile,labelname = getLabelFile(label,sub,hemi)
            labelfiles.append(labelfile)
            labelnames.append(labelname)
        annot_file = label2Annot(sub,hemi,labelnames,labelfiles,annot_name)
        annot2label(sub,hemi,annot_file)

        label_files = glob.glob('%s/%s/label/sulcal_prob_mpm/%s.*.label'%(os.environ['SUBJECTS_DIR'],sub,hemi))
        checkForOverlap(label_files)
