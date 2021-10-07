if isfield(output,'fname')
    % replace unallowed characters
    sulcnames = output.list_sulc;
    sulcnames = strrep(sulcnames,'&','');
    sulcnames = strrep(sulcnames,'-','');
    sulcnames = [strcat('lh_',sulcnames)s]; %% change to lh
    
    % write sulci width csv
    tbl = array2table(output.sulci_width,'VariableNames',sulcnames);
    tbl = [table(output.list_subject','VariableNames',{'subID_calcSulc'}) tbl];
    writetable(tbl,[output.fname '_sulciwidth_sulc_hemi.csv']); %% change sulc and hemi according to list in HCP_wrapper
    
    % write sulci depth csv
    tbl = array2table(output.sulci_depth,'VariableNames',sulcnames);
    tbl = [table(output.list_subject','VariableNames',{'subID_calcSulc'}) tbl];
    writetable(tbl,[output.fname '_sulcidepth_insula.csv']); %% change sulc and hemi according to list in HCP_wrapper
end