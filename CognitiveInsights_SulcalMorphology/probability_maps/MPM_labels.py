print('Begin MPM python script')
# imports
import numpy as np
import nibabel as nib
import nilearn
from nilearn import image, plotting, surface
import scipy
import matplotlib.pyplot as plt
import os, sys
import glob
import seaborn as sns
import pandas as pd





left_out_sub = sys.argv[1]
subject_dir = sys.argv[2]
fsaverage_space_labels = sys.argv[3] 



# load subjects from list of processed subs
subjects = np.genfromtxt(sys.argv[4], dtype = str)



# Project identifier for the MPM file
project_id = sys.argv[5]


#subjects = np.delete(subjects, left_out_sub_index)

# Uncomment to review inputs from .sh 
#print('ARGV 1: ', sys.argv[1])
#print('ARGV 2: ', sys.argv[2])
#print('ARGV 3: ', sys.argv[3])
#print('ARGV 4: ', sys.argv[4])
#print('ARGV 5: ', sys.argv[5])
#print('ARGV 6: ', sys.argv[6])

prediction_fsaverage_space_labels = sys.argv[6]

subjects = np.genfromtxt(sys.argv[7], dtype = str)

middle_frontal_label_names = np.array(sys.argv[8:])

print('Assembling MPMs for, ', middle_frontal_label_names)
print('included subjects: {}\n\n\n'.format(subjects))


# define useful functions
def read_label(label_name):
    """
    Reads a freesurfer-style .label file (5 columns)

    Parameters
    ----------
    label_name: str

    Returns
    -------
    vertices: index of the vertex in the label np.array [n_vertices]
    RAS_coords: columns are the X,Y,Z RAS coords associated with vertex number in the label, np.array [n_vertices, 3]

    """

    # read label file, excluding first two lines of descriptor
    df_label = pd.read_csv(label_name,skiprows=[0,1],header=None,names=['vertex','x_ras','y_ras','z_ras','stat'],delimiter='\s+')

    vertices = np.array(df_label.vertex)
    RAS_coords = np.empty(shape = (vertices.shape[0], 3))
    RAS_coords[:,0] = df_label.x_ras
    RAS_coords[:,1] = df_label.y_ras
    RAS_coords[:,2] = df_label.z_ras

    return vertices, RAS_coords
def read_label_stat(label_name):
    """
    Reads a freesurfer-style .label file (5 columns)

    Parameters
    ----------
    label_name: str

    Returns
    -------
    vertices: index of the vertex in the label np.array [n_vertices]
    RAS_coords: columns are the X,Y,Z RAS coords associated with vertex number in the label, np.array [n_vertices, 3]

    """

    # read label file, excluding first two lines of descriptor
    df_label = pd.read_csv(label_name,skiprows=[0,1],header=None,names=['vertex','x_ras','y_ras','z_ras','stat'],delimiter='\s+')

    vertices = np.array(df_label.vertex)
    RAS_coords = np.empty(shape = (vertices.shape[0], 3))
    RAS_coords[:,0] = df_label.x_ras
    RAS_coords[:,1] = df_label.y_ras
    RAS_coords[:,2] = df_label.z_ras
    stat = df_label.stat

    return vertices, RAS_coords, stat



### make probability maps
# loop through labels, and for each count the % of subjects at each vertex

vertices_dict_lh = {}
vertices_dict_rh = {}

print('Making probability maps for ', left_out_sub, '\n\n')

for i, middle_frontal_label in enumerate(middle_frontal_label_names):

    label_vertices_lh = np.empty(shape=0,dtype=int)
    label_RAS_lh = np.empty(shape=(0,3),dtype=int)
    label_vertices_rh = np.empty(shape=0,dtype=int)
    label_RAS_rh = np.empty(shape=(0,3),dtype=int)
    

    # Make directory for held-out subject, save out probability maps from all other subjects there

    prob_label_path_lh = fsaverage_space_labels + '/prob_maps/{}/lh.{}_PROB_{}.label'.format(left_out_sub, project_id, middle_frontal_label)
    prob_label_path_rh = fsaverage_space_labels + '/prob_maps/{}/rh.{}_PROB_{}.label'.format(left_out_sub, project_id,  middle_frontal_label)
    


    
    ## Create left hemisphere _PROB_ labels
        
    for sub in subjects:
        try:
        
            label_path_lh = prediction_fsaverage_space_labels + 'projected_labels/{}/{}.lh.{}.label'.format(middle_frontal_label, sub, middle_frontal_label)
            
            vertices_lh, RAS_lh = read_label(label_path_lh)
            label_vertices_lh = np.append(label_vertices_lh,vertices_lh)
            label_RAS_lh = np.append(label_RAS_lh,RAS_lh,axis=0)
            
            unique_vertices,indices_vertices,counts_vertices=np.unique(label_vertices_lh,return_index=True,return_counts=True)
            # index only the RAS coords for unique vertices
            unique_RAS = label_RAS_lh[indices_vertices,:]
            # get probabilities at each vertex
            
            prob_vertices = (counts_vertices)/(subjects.shape[0])
            

            # make probabilistic label array for label file
            prob_array = np.zeros(shape=(unique_vertices.shape[0],5),dtype=float)
            prob_array[:,0] = unique_vertices
            prob_array[:,1:4] = unique_RAS
            prob_array[:,-1] = prob_vertices

            np.savetxt(prob_label_path_lh, prob_array, fmt='%-2d  %2.3f  %2.3f  %2.3f %1.10f')

            # edit first two lines of label file to match Freesurfer
            
            f = open(prob_label_path_lh, 'r')
            edit_f = f.read()
            f.close()
            f = open(prob_label_path_lh, 'w')
            f.write('#!ascii label  , from subject fsaverage vox2ras=TkReg\n{}\n'.format(unique_vertices.shape[0]))
            f.write(edit_f)
            f.close()
            

        except Exception:
            pass
            #print('Left hemisphere ', label_path_lh, ' does not exist for ', sub)

    print('Left Hemishere ', middle_frontal_label, 'prob_map written for ', left_out_sub, ' \n\n\n') 

    ## Create right hemisphere _PROB_ labels
    
    for sub in subjects:
        try:
            label_path_rh = prediction_fsaverage_space_labels + 'projected_labels/{}/{}.rh.{}.label'.format(middle_frontal_label, sub ,middle_frontal_label)
            
            vertices_rh, RAS_rh = read_label(label_path_rh)
            label_vertices_rh = np.append(label_vertices_rh,vertices_rh)
            label_RAS_rh = np.append(label_RAS_rh,RAS_rh,axis=0)

            # repeat process for RH
            # get unique vertices and their index and counts
            unique_vertices,indices_vertices,counts_vertices=np.unique(label_vertices_rh,return_index=True,return_counts=True)
            # index only the RAS x`coords for unique vertices
            unique_RAS = label_RAS_rh[indices_vertices,:]
            # get probabilities at each vertex
            prob_vertices = (counts_vertices)/(subjects.shape[0])
            

            # make probabilistic label array for label file
            prob_array = np.zeros(shape=(unique_vertices.shape[0],5),dtype=float)
            prob_array[:,0] = unique_vertices
            prob_array[:,1:4] = unique_RAS
            prob_array[:,-1] = prob_vertices
            np.savetxt(prob_label_path_rh, prob_array, fmt='%-2d  %2.3f  %2.3f  %2.3f %1.10f')
            
            # edit first two lines of label file to match Freesurfer
            f = open(prob_label_path_rh, 'r')
            edit_f = f.read()
            f.close()
            f = open(prob_label_path_rh, 'w')
            f.write('#!ascii label  , from subject fsaverage vox2ras=TkReg\n{}\n'.format(unique_vertices.shape[0]))
            f.write(edit_f)
            f.close()
        except Exception:
            pass
            

    print('Right Hemishere ', middle_frontal_label, 'prob_map written for ', left_out_sub, ' \n\n\n') 

    
# make max probability maps
# loop through each vertex in the hemisphere and determine which labels contain the vertex
# for each vertex, see which labels contain the vertex, and then find the label with the highest stat at that vertex,
# and add that vertex number and stat to the new MPM version of the label

# ### left hemi
# load probability maps and enter into dict
prob_maps_vertices = {}
prob_maps_RAS = {}
prob_maps_stat = {}

MPM_vertices = {}
MPM_RAS = {}
MPM_stat = {}

print('Creating Left Hemisphere MPMs \n\n\n')

# loop through labels, load prob map and make empty values in MPM dicts
for i, middle_frontal_label in enumerate(middle_frontal_label_names):
    try:
        #load the prob mpa for the given label
        vertices, RAS, stat = read_label_stat(fsaverage_space_labels + 'prob_maps/{}/lh.{}_PROB_{}.label'.format(left_out_sub, project_id, middle_frontal_label))
        prob_maps_vertices[middle_frontal_label] = vertices
        prob_maps_RAS[middle_frontal_label] = RAS
        prob_maps_stat[middle_frontal_label] = stat

        MPM_vertices[middle_frontal_label] = np.empty(0)
        MPM_RAS[middle_frontal_label] =  np.empty(shape=(0,3))
        MPM_stat[middle_frontal_label] = np.empty(0)
    except Exception:
        pass
        #load lh cortex vertices
vertices_lh, RAS_lh = read_label('/home/weiner/data/fsaverage/label/lh.cortex.label')
    

# loop through lh cortex vertices
for vtx in vertices_lh:

    labels_with_vtx = np.empty(0)
    vertices_with_vtx = np.empty(0)
    RAS_with_vtx = np.empty(shape=(0,3))
    stat_with_vtx = np.empty(0)

    for label_prob, vertices_prob in prob_maps_vertices.items():
        # if vtx from cortex is in vertices of prob map, add name to list of labels with vtx, add stat value
        match_idx = np.where(vertices_prob == vtx)
        # if vertex is in probability map
        if match_idx[0].shape[0] > 0:
            # get vertex, RAS, and stat values for the given vertex that is in prob map
            vertices_prob_idx = match_idx[0][0]
            labels_with_vtx = np.append(labels_with_vtx, label_prob)
            vertices_with_vtx = np.append(vertices_with_vtx, vertices_prob[vertices_prob_idx])
            RAS_with_vtx = np.concatenate((RAS_with_vtx, np.reshape(prob_maps_RAS[label_prob][vertices_prob_idx,:],(1,3))),axis=0)
            stat_with_vtx = np.append(stat_with_vtx, prob_maps_stat[label_prob][vertices_prob_idx])

     # if vertex was in at least one prob map, get the max value and add to MPM file
        # also required to have a probability of 0.33 or higher
    if (labels_with_vtx.shape[0] > 0):

        if (np.max(stat_with_vtx) > 1/3):

            max_idx = np.argmax(stat_with_vtx)

            label_max = labels_with_vtx[max_idx]
            RAS_max = RAS_with_vtx[max_idx,:]
            stat_max = stat_with_vtx[max_idx]

            MPM_vertices[label_max] = np.append(MPM_vertices[label_max], vtx)
            MPM_RAS[label_max] = np.concatenate((MPM_RAS[label_max], np.reshape(RAS_max,(1,3))),axis=0)
            MPM_stat[label_max] = np.append(MPM_stat[label_max], stat_max)

# save out each entry MPM as a separate label file in Freesurfer

for i, middle_frontal_label in enumerate(middle_frontal_label_names):

    MPM_path = fsaverage_space_labels + 'prob_maps/{}/lh.{}_PROB_MPM_{}.label'.format(left_out_sub, project_id, middle_frontal_label)
    try:
        # make probabilistic label array for albel file
        prob_array = np.zeros(shape=(MPM_vertices[middle_frontal_label].shape[0],5),dtype=float)
        prob_array[:,0] = MPM_vertices[middle_frontal_label]
        prob_array[:,1:4] = MPM_RAS[middle_frontal_label]
        prob_array[:,-1] = MPM_stat[middle_frontal_label]
        np.savetxt(MPM_path, prob_array, fmt='%-2d  %2.3f  %2.3f  %2.3f %1.10f')

        # edit first two lines of label file to match Freesurfer
        f = open(MPM_path, 'r')
        edit_f = f.read()
        f.close()
        f = open(MPM_path, 'w')
        f.write('#!ascii label  , from subject fsaverage vox2ras=TkReg\n{}\n'.format(MPM_vertices[middle_frontal_label].shape[0]))
        f.write(edit_f)
        f.close()
    except Exception:
        pass
print('Left Hemisphere PROB MPMs written for ', left_out_sub,' \n\n\n') 

# save out each entry MPM as a separate binary label file in Freesurfer

for i, middle_frontal_label in enumerate(middle_frontal_label_names):

    MPM_path = fsaverage_space_labels + 'prob_maps/{}/lh.{}_PROB_MPM_binary_{}.label'.format(left_out_sub, project_id ,middle_frontal_label)
    try:
        # make probabilistic label array for albel file
        prob_array = np.zeros(shape=(MPM_vertices[middle_frontal_label].shape[0],5),dtype=float)
        prob_array[:,0] = MPM_vertices[middle_frontal_label]
        prob_array[:,1:4] = MPM_RAS[middle_frontal_label]
        prob_array[:,-1] = 1
        np.savetxt(MPM_path, prob_array, fmt='%-2d  %2.3f  %2.3f  %2.3f %1.10f')

        # edit first two lines of label file to match Freesurfer
        f = open(MPM_path, 'r')
        edit_f = f.read()
        f.close()
        f = open(MPM_path, 'w')
        f.write('#!ascii label  , from subject fsaverage vox2ras=TkReg\n{}\n'.format(MPM_vertices[middle_frontal_label].shape[0]))
        f.write(edit_f)
        f.close()

    except Exception:
        pass
    
print('Left Hemisphere PROB Binary MPMs written for ', left_out_sub, '\n\n\n')

# ### right hemi
# load probability maps and enter into dict
prob_maps_vertices = {}
prob_maps_RAS = {}
prob_maps_stat = {}

MPM_vertices = {}
MPM_RAS = {}
MPM_stat = {}

print('Creating Right Hemisphere MPMs \n\n\n')
# loop through labels, load prob map and make empty values in MPM dicts
for i, middle_frontal_label in enumerate(middle_frontal_label_names):
    #load the prob mpa for the given label
    try:
        vertices, RAS, stat = read_label_stat(fsaverage_space_labels + 'prob_maps/{}/rh.{}_PROB_{}.label'.format(left_out_sub, project_id, middle_frontal_label))
        prob_maps_vertices[middle_frontal_label] = vertices
        prob_maps_RAS[middle_frontal_label] = RAS
        prob_maps_stat[middle_frontal_label] = stat

        MPM_vertices[middle_frontal_label] = np.empty(0)
        MPM_RAS[middle_frontal_label] =  np.empty(shape=(0,3))
        MPM_stat[middle_frontal_label] = np.empty(0)
    except Exception:
        pass

#load rh cortex vertices
vertices_rh, RAS_rh = read_label('/home/weiner/data/fsaverage/label/rh.cortex.label')

# loop through rh cortex vertices
for vtx in vertices_rh:

    labels_with_vtx = np.empty(0)
    vertices_with_vtx = np.empty(0)
    RAS_with_vtx = np.empty(shape=(0,3))
    stat_with_vtx = np.empty(0)

    for label_prob, vertices_prob in prob_maps_vertices.items():
        # if vtx from cortex is in vertices of prob map, add name to list of labels with vtx, add stat value
        match_idx = np.where(vertices_prob == vtx)
        # if vertex is in probability map
        if match_idx[0].shape[0] > 0:
            # get vertex, RAS, and stat values for the given vertex that is in prob map
            vertices_prob_idx = match_idx[0][0]
            labels_with_vtx = np.append(labels_with_vtx, label_prob)
            vertices_with_vtx = np.append(vertices_with_vtx, vertices_prob[vertices_prob_idx])
            RAS_with_vtx = np.concatenate((RAS_with_vtx, np.reshape(prob_maps_RAS[label_prob][vertices_prob_idx,:],(1,3))),axis=0)
            stat_with_vtx = np.append(stat_with_vtx, prob_maps_stat[label_prob][vertices_prob_idx])

     # if vertex was in at least one prob map, get the max value and add to MPM file
        # also required to have a probability of 0.33 or higher
    if (labels_with_vtx.shape[0] > 0):

        if (np.max(stat_with_vtx) > 1/3):

            max_idx = np.argmax(stat_with_vtx)

            label_max = labels_with_vtx[max_idx]
            RAS_max = RAS_with_vtx[max_idx,:]
            stat_max = stat_with_vtx[max_idx]

            MPM_vertices[label_max] = np.append(MPM_vertices[label_max], vtx)
            MPM_RAS[label_max] = np.concatenate((MPM_RAS[label_max], np.reshape(RAS_max,(1,3))),axis=0)
            MPM_stat[label_max] = np.append(MPM_stat[label_max], stat_max)

# save out each entry MPM as a separate label file in Freesurfer

for i, middle_frontal_label in enumerate(middle_frontal_label_names):

    MPM_path = fsaverage_space_labels + 'prob_maps/{}/rh.{}_PROB_MPM_{}.label'.format(left_out_sub, project_id, middle_frontal_label)
    try:
        # make probabilistic label array for albel file
        prob_array = np.zeros(shape=(MPM_vertices[middle_frontal_label].shape[0],5),dtype=float)
        prob_array[:,0] = MPM_vertices[middle_frontal_label]
        prob_array[:,1:4] = MPM_RAS[middle_frontal_label]
        prob_array[:,-1] = MPM_stat[middle_frontal_label]

        np.savetxt(MPM_path, prob_array, fmt='%-2d  %2.3f  %2.3f  %2.3f %1.10f')

        # edit first two lines of label file to match Freesurfer
        f = open(MPM_path, 'r')
        edit_f = f.read()
        f.close()
        f = open(MPM_path, 'w')
        f.write('#!ascii label  , from subject fsaverage vox2ras=TkReg\n{}\n'.format(MPM_vertices[middle_frontal_label].shape[0]))
        f.write(edit_f)
        f.close()
    except Exception:
        pass

print('Right Hemisphere PROB MPMs written for ', left_out_sub, '\n\n\n')


# save out each entry MPM as a separate binary label file in Freesurfer

for i, middle_frontal_label in enumerate(middle_frontal_label_names):

    MPM_path = fsaverage_space_labels + 'prob_maps/{}/rh.{}_PROB_MPM_binary_{}.label'.format(left_out_sub, project_id, middle_frontal_label)
    try:
        # make probabilistic label array for label file
        prob_array = np.zeros(shape=(MPM_vertices[middle_frontal_label].shape[0],5),dtype=float)
        prob_array[:,0] = MPM_vertices[middle_frontal_label]
        prob_array[:,1:4] = MPM_RAS[middle_frontal_label]
        prob_array[:,-1] = 1
        np.savetxt(MPM_path, prob_array, fmt='%-2d  %2.3f  %2.3f  %2.3f %1.10f')

        # edit first two lines of label file to match Freesurfer
        f = open(MPM_path, 'r')
        edit_f = f.read()
        f.close()
        f = open(MPM_path, 'w')
        f.write('#!ascii label  , from subject fsaverage vox2ras=TkReg\n{}\n'.format(MPM_vertices[middle_frontal_label].shape[0]))
        f.write(edit_f)
        f.close()
    except Exception:
        pass
print('Right Hemisphere PROB Binary MPMs written for ', left_out_sub, '\n\n\n')

