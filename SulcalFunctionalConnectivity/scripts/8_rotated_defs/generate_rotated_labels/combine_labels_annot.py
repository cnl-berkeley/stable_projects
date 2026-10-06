import os,sys
import numpy as np
import seaborn as sns

def getLabelFile(label,sub,hemi):
    labelfile_found = ''
    labelname_found = ''
    if type(label) is list: # pick the first label in list that exits
        label_found = False
        for lab in label:
            labelfile = '%s/%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,lab)
            if os.path.exists(labelfile) and not label_found:
                labelfile_found = labelfile
                label_found = True
                labelname_found = label[0]
    else: # look for specific label
        labelfile = '%s/%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,label)
        if os.path.exists(labelfile):
            labelfile_found = labelfile
            labelname_found = label
 
    print(label,labelfile_found,labelname_found)
    return labelfile_found,labelname_found


def getLUT(labels):
    colors = sns.color_palette("hls", len(labels))
    colors = colors
    lut = '%s/label_LUT.txt'%os.getcwd()
    of = open(lut,'w')
    of.write('#No. Label Name:                            R   G   B   A\n')
    of.write('0   Unknown                                 0   0   0   0\n')
    for l in np.arange(len(labels)):
        of.write('%d   %s                                 %d   %d   %d   0\n'%(l+1,labels[l],colors[l][0]*255,colors[l][1]*255,colors[l][2]*255))
    of.close()

    return lut


def label2Annot(sub,hemi,labelnames,label_files):
    from nipype.interfaces.freesurfer import Label2Annot
    l2a = Label2Annot()
    l2a.inputs.hemisphere = hemi
    l2a.inputs.subject_id = sub
    l2a.inputs.in_labels = label_files
    l2a.inputs.orig = '%s/%s/surf/%s.pial'%(os.environ['SUBJECTS_DIR'],sub,hemi)
    l2a.inputs.out_annot = 'LPFC+LPC' #.annot'
    l2a.inputs.color_table = getLUT(labelnames)
    print(l2a.cmdline)
    l2a.run()


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
    else:
        lpfc_labels = ['ifs',['painfs_any','painfs_combined'],'pmfs_a','pmfs_i','pmfs_p','sfs_a','sfs_p'] + ['prts','lfms','aalf']
        lpar_labels = [['slos1','slocs-v','SLOS'],'sB','pips','mTOS',['iTOS','ITOS','lTOS'],'IPS-PO','IPS','cSTS1','cSTS2','cSTS3','aipsJ']
        mpar_labels = []
        labels = lpfc_labels + lpar_labels + mpar_labels

    if sub != 'fsaverage_prob':
        sub = 'sub-%s'%sub

    for hemi in ['lh','rh']:
        labelfiles = []
        labelnames = []
        for label in labels:
            labelfile,labelname = getLabelFile(label,sub,hemi)
            labelfiles.append(labelfile)
            labelnames.append(labelname)
        label2Annot(sub,hemi,labelnames,labelfiles)

