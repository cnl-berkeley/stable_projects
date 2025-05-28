import os,sys
from numpy import *

def resampleLabel(labelname,sub,hemi,outdir):
    # resample label to 32k space

    # pick the same surfaces used for the fc matrix (or be sure to use the same surface e.g. midthickness!)
    fmriprep_derivatives_dir = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/bids/derivatives_v7/sub-%s/ses-%s'%(sub,sub[-1])
    sdir='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/cortex/fc_32k232k_0mm-bp-relmatch/sub-%s/surf_32k/gifti'%(sub,sub)
    highres_surface = '%s/anat/sub-%s_ses-%s_hemi-%s_midthickness.surf.gii'%(fmriprep_derivatives_dir,sub,sub[-1],hemi[0].upper())
    lowres_surface = '%s/%s-32k.surf.gii'%(sdir,os.path.basename(highres_surface).split('.')[0])
    highres_sphere = '%s/%s.sphere.reg.surf.gii'%(sdir,hemi)
    lowres_sphere = '%s/sphere.32k.%s.surf.gii'%(sdir,hemi.upper()[0])

    # convert label to gii
    highres_surface_bin = '%s/sub-%s/surf/%s.white'%(os.environ['SUBJECTS_DIR'],sub,hemi)
    highres_label_bin = '%s/sub-%s/label/%s.%s.label'%(os.environ['SUBJECTS_DIR'],sub,hemi,labelname)
    highres_label = '%s/%s.%s.label.gii'%(outdir,labelname,hemi[0].upper())
    str = 'mris_convert --label %s %s %s %s'%(highres_label_bin,labelname,highres_surface_bin,highres_label)
    os.system(str)

    # resample
    lowres_label = '%s/%s-32k.%s.label.gii'%(outdir,labelname,hemi[0].upper())
    import nipype.interfaces.workbench.base as wb
    str = 'wb_command -label-resample %s %s %s ADAP_BARY_AREA -area-surfs %s %s %s'%(highres_label,highres_sphere,lowres_sphere,highres_surface,lowres_surface,lowres_label)
    cmd = wb.WBCommand(str)
    cmd.run()

    sdir='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/cortex/fc_32k232k_0mm-bp-relmatch/sub-%s/surf_32k/fs'%(sub,sub)
    lowres_surface_bin = '%s/%s-32k'%(sdir,os.path.basename(highres_surface).split('.')[0])

    return lowres_surface_bin,lowres_label


def getAdjMatrix(lowres_surface_bin,lowres_label,outdir):
    import nibabel as nib
    # calculate adjacency matrix between 32k surface nodes within the given label
    # NOTE: medial wall is not explicitly masked from being included in shortest paths

    # load surface, using inflated as in the nilearn example
    # https://nilearn.github.io/dev/auto_examples/02_decoding/plot_haxby_searchlight_surface.html#surface-bold-response
    vertices, triangles = nib.freesurfer.io.read_geometry(lowres_surface_bin)
    # https://stackoverflow.com/questions/66757983/big-endian-is-not-supported-in-little-endian-compiler-python-how-to-convert-to
    triangles = triangles.byteswap().newbyteorder('=')

    # load indices of label on surface
    giidata = nib.load(lowres_label)
    label_vec = squeeze(asarray([x.data for x in giidata.darrays]))
    label_inds = where(label_vec > 0)[0]
    label_inds = label_inds.astype(int32)
    print('label_inds',len(label_inds))

    # calculate pairwise geodesic distances
    import gdist
    print('calculating distance matrix')
    adj = gdist.distance_matrix_of_selected_points(vertices,triangles,label_inds,50.0) # 5cm

    adj_rounded = around(adj,decimals=1)
    adj_rounded = adj_rounded*10
    adj_uint = adj_rounded.astype(uint16)

    print('saving distance matrix 2')
    adj_full = adj_uint.toarray()
    adj = adj_full[ix_(label_inds,label_inds)]
    print(amax(adj))
    print('%s/adj3_uint16-%s-%s.npz'%(outdir,os.path.basename(lowres_label).split('.')[0],hemi))
    savez_compressed('%s/adj3_uint16-%s-%s.npz'%(outdir,os.path.basename(lowres_label).split('.')[0],hemi),adj)


if __name__ == '__main__':

    os.environ['SUBJECTS_DIR'] = '/home/weiner/Nora_PFCSulci/Projects/NORA_relmatch_funcNeuroAnat/data/subjects_v7'

    sub = sys.argv[1]
    hemi = sys.argv[2]
    labelname = 'cortex'

    outdir='/home/weiner/shakki/NotBackedUp/spatial_autocorrelation/%s/cortex'%sub
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    print(outdir)

    lowres_surface_bin,lowres_label = resampleLabel(labelname,sub,hemi,outdir)
    print('*')
    getAdjMatrix(lowres_surface_bin,lowres_label,outdir)

