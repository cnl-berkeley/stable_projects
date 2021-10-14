import pandas as pd
import numpy as np
import os
import sys


def get_average_vals(subjects_dir, subject_list, label_name, project_dir):
    ## Find average RAS coord
    ## SAVES TO DF in SUBJECT DIR under name label_name_posterior_axis_vals.csv

    ## RUN FROM COMMAND LINE AS:
    ## $ python RAS_means.py $subjects_dir $subjects_list $label_name $project_dir
    

    subs = []
    label_df = pd.DataFrame(columns=['sub', 'hemi'])
    centroid = ['centroid']
    hemis = ['lh', 'rh']
    RAS_columns = ['R', 'A', 'S']

    with open(subject_list) as sample:
        for line in sample:
            subs.append(line.rstrip())

    sub_hemi = [[a, b] for a in hemis for b in subs]
    
    for i, sub_hemi_pair in enumerate(sub_hemi):
        for num, col in enumerate(RAS_columns):
            hemi = sub_hemi_pair[0]
            sub = sub_hemi_pair[1]

            label_df.set_value(i, 'sub', sub)
            label_df.set_value(i, 'hemi', hemi)
            

            label_path = os.path.join(os.path.sep, subjects_dir, sub, 'label/{}.{}.label'.format(hemi, label_name))
            print(label_path)
            try:
                label_axis = np.loadtxt(label_path, skiprows=2)[:,num+1]
        
                print('Label found ...')
                label_mean = np.mean(label_axis)
                

                label_df.set_value(i, '{}_{}_coord'.format(label_name, col), label_mean)
                print('Finished ', sub, hemi, col)
            except Exception:
                label_df.set_value(i, '{}_{}_coord'.format(label_name, col), np.nan)
                print('Finished Nan', sub, hemi, col)
                
                 
    
    label_df.dropna()
    df_filepath = os.path.join(os.path.sep, project_dir, 'Network_{}_mean_axis_vals.csv'.format(label_name))
    label_df.to_csv(df_filepath)
    print('DataFrame written at', df_filepath)


if __name__=='__main__':
    args = [str(i) for i in sys.argv]
        
    get_average_vals(args[1], args[2], args[3], args[4])
      
