function depth = calcSulc_depth(options,subject_hemi,mesh)

% identify fundus (valley) for depth
try
    % find vertices for lowest points in map, using sulcal map
    % lowest 100 vertices
    [v,i]=sort(subject_hemi.sulcmap(mesh.label_v));
    if length(mesh.label_v(i)) > 100
        fundus_v = mesh.label_v(i((end-99):end));
      
    else
        fundus_v = mesh.label_v(i((end-10):end));
    
    end
    
    % calculate distance between vertex on pial surface and gyrif surface vertex
    for idx = 1:length(fundus_v)
        v_xyz = subject_hemi.pial_v(fundus_v(idx),:);
        depth(idx) = min(pdist2(v_xyz,subject_hemi.gyrif_v));
    end
    
    depth = median(depth);
catch
    depth = NaN;
disp(depth); 
end