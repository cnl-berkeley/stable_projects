function new_cius = ciu_value_fix(ciu_file)
% clean up fancy

cius = load(ciu_file);

new_cius = zeros(size(cius));
new_cius(:,end) = cius(:,end);

for i = 1:(length(cius(2,:))-1)
    j = length(cius(1,:))-i;
    
    e = 1;
    uvals = unique(cius(:,j));
    uvals_prev = unique(new_cius(:,j+1));
    
    overlap_order = zeros(length(uvals),length(uvals_prev));
    % compare a cluster against all clusters of previous round
    % get "correspondence" matrix
    for u = 1:length(uvals)
        cluster_inds = find(cius(:,j) == uvals(u)); %indices of nodes in this cluster
        for v = 1:length(uvals_prev)
            inds1 = find(new_cius(cluster_inds,j+1) == uvals_prev(v));
            overlap_order(u,v) = length(inds1);
        end
    end
    
    proc_uvals = zeros(length(uvals),1);
    set_uvals_prev = zeros(length(uvals_prev),1);
    
    % reorder, start from most overlapping clusters
    for u=1:length(uvals)
        [i1,i2] = find(overlap_order == max(max(overlap_order))); % find cluster associated with most overlapping nodes
        if max(max(overlap_order)) == 0 % no overlap defined (cleanly split cluster, smaller part)
            if proc_uvals(u) == 0
                inds = find(cius(:,j) == uvals(i1(1))); % find indices of the cluster
                
                % determine new value, from previous iterations if available
                new_val = -1;
                n = j+1;
                while n < length(cius(1,:))
                    val = mode(new_cius(inds,n));
                    if length(find(unique(new_cius(:,j)) == val)) == 0 % value not used this iteration
                        new_val = val;
                    end
                    n = n+1;
                end
                if new_val > 0
                    new_cius(inds,j) = new_val;
                else % not possible
                    new_cius(inds,j) = max([uvals;uvals_prev])+e; %closest_val+1;
                    e = e+1;
                end
                overlap_order(i1,:) = 0; % not on option for other clusters
            end
            
        else
            t=1;
            while t <= length(i1) % loop cases with the same max overlap
                inds = find(cius(:,j) == uvals(i1(t))); % find indices of the cluster
                % set the new value
                if max(max(overlap_order)) > 0
                    if set_uvals_prev(i2(t)) == 0 % catches if a new cluster is formed
                        new_cius(inds,j) = uvals_prev(i2(t));
                        set_uvals_prev(i2(t)) = 1;
                    else
                        % determine new value, from previous iterations if available
                        new_val = -1;
                        n = j+1;
                        while n < length(cius(1,:))
                            val = mode(new_cius(inds,n));
                            if length(find(unique(new_cius(:,j)) == val)) == 0 % value not used this iteration
                                new_val = val;
                            end
                            n = n+1;
                        end
                        if new_val > 0
                            new_cius(inds,j) = new_val;
                        else % not possible
                            new_cius(inds,j) = max([uvals;uvals_prev])+e; %closest_val+1;
                            e = e+1;
                        end
                        %new_cius(inds,j) = max([uvals;uvals_prev])+e;
                        %e = e+1;
                    end
                else
                    if proc_uvals(u) == 0
                        % determine new value, from previous iterations if available
                        new_val = -1;
                        n = j+1;
                        while n < length(cius(1,:))
                            val = mode(new_cius(inds,n));
                            if length(find(unique(new_cius(:,j)) == val)) == 0 % value not used this iteration
                                new_val = val;
                            end
                            n = n+1;
                        end
                        if new_val > 0
                            new_cius(inds,j) = new_val;
                        else % not possible
                            new_cius(inds,j) = max([uvals;uvals_prev])+e; %closest_val+1;
                            e = e+1;
                        end
                    end
                end
                
                % remove the processed
                overlap_order(i1(t),:) = 0; % not on option for other clusters
                proc_uvals(i1(t)) = 1;
                
                t = t+1;
            end
        end
    end

    [path,name,ext] = fileparts(ciu_file);
    csvwrite(sprintf('%s/%s-fixed.txt',path,name),new_cius);
end
