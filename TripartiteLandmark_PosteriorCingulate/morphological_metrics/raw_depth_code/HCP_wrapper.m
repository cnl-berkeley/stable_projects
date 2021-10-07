

subject_dir = '/home/weiner/HCP/subjects/';


%% STEPS: 
% Update subject_dir
% Update lh_sub_categories
% be sure field names match .txt for each subject category
% Update sulcal names to match
% Change calcsave hemisphere on commented line 6


lh_sub_categories = struct();
lh_sub_categories.lh_2 = {'MCGS', 'sbps', 'prcus1', 'prcus2', 'prcus3', 'prculs', 'POS', '2'};
lh_sub_categories.lh_12 = {'MCGS', 'sbps', 'prcus1', 'prcus2', 'prcus3', 'prculs', 'POS', '1', '2' };
lh_sub_categories.lh_23 = {'MCGS', 'sbps', 'prcus1', 'prcus2', 'prcus3', 'prculs', 'POS', '2', '3'};
lh_sub_categories.lh_123 = {'MCGS', 'sbps', 'prcus1', 'prcus2', 'prcus3', 'prculs', 'POS', '1', '2', '3'};
lh_sub_categories.lh_1234 = {'MCGS', 'sbps', 'prcus1', 'prcus2', 'prcus3', 'prculs', 'POS', '1', '2', '3', '4'};
%lh_sub_categories.lh_sps = {'sps'};

rh_sub_categories = struct();
%rh_sub_categories.rh_2 = {'MCGS', 'sbps', 'prcus1', 'prcus2', 'prcus3', 'prculs', 'POS', '2',};
%rh_sub_categories.rh_12 = {'MCGS', 'sbps', 'prcus1', 'prcus2', 'prcus3', 'prculs', 'POS', '1', '2' };
%%rh_sub_categories.rh_23 = {'MCGS', 'sbps', 'prcus1', 'prcus2', 'prcus3', 'prculs', 'POS', '2', '3'};
%rh_sub_categories.rh_123 = {'MCGS', 'sbps', 'prcus1', 'prcus2', 'prcus3', 'prculs', 'POS', '1', '2', '3'};
%rh_sub_categories.rh_234 = {'MCGS', 'sbps', 'prcus1', 'prcus2', 'prcus3', 'prculs', 'POS', '2', '3', '4'};
%rh_sub_categories.rh_12345 = {'sbps', 'prcus1', 'prcus2', 'prcus3', 'prculs', 'POS', '1', '2', '3', '4', '5'};
rh_sub_categories.rh_sps = {'sps'};



lh_fieldn = fieldnames(lh_sub_categories);
rh_fieldn = fieldnames(rh_sub_categories);

%% Uncomment all of left hemisphere
% ALSO change lh / rh in calcSulc.m with whichever you're running

%% Loop through each of the sulcal lists
for i=1:numel(fieldnames(lh_sub_categories))
    
    subject_group = lh_fieldn(i);
    subject_group_str = lh_fieldn{i};
    
    %%% UPDAT FILE PATH
    file_path = strcat(fullfile('/home/weiner/HCP/projects/ifrms_HCP/divided_subs/loc_lists',subject_group),'.txt');
    
 
    subjects_file = fopen(file_path{1}, 'r');
   
    subjects_cell = textscan(subjects_file,'%s','delimiter','\n');
    
    subjects = subjects_cell{1}(1:end);
    %% loop through each of the subjects in each sulcal category
   for n=1:length(subjects)
        subject = {subjects{n}};
        disp(file_path)
        disp(subject)
        options.list_sulc = getfield(lh_sub_categories, subject_group_str);
        
        options.estimateWidth = 1;
        options.estimateDepth = 1;
        options.hemis = 1;
        
        output = calcSulc(subject,subject_dir,options);
        output.fname = sprintf('%s_calcsulc_values',subject{1});
        
        %% SAVED TO CURRENT WORKING DIRECTORY
        save(sprintf('%s_lh_calcsulc_values_v2.mat',subject{1}),'output');
        
        %%% CSV SAVED IN calvcSulc_save.m
        calcSulc_save
   end
end

%disp('End left hemisphere')


%for i=1:numel(fieldnames(rh_sub_categories))
    
    
    %subject_group = rh_fieldn(i);
    %subject_group_str = rh_fieldn{i};
    
    %%% UPDATE FILE PATH
    %file_path = strcat(fullfile('/home/weiner/HCP/projects/ifrm_HCP/divided_subs/loc_lists',subject_group),'.txt');
    %disp(file_path)
    %subjects_file = fopen(file_path{1}, 'r');
   
    %subjects_cell = textscan(subjects_file,'%s','delimiter','\n');
    
    %subjects = subjects_cell{1}(1:end);
    
   %for n=1:length(subjects)
        %subject = {subjects{n}};
        %options.list_sulc = getfield(rh_sub_categories, subject_group_str);
        %options.estimateWidth = 1;
        %options.estimateDepth = 1;
        %options.hemis = 0;
        %disp(subject)
        
        %output = calcSulc(subject,subject_dir,options);
        %output.fname = sprintf('%s_calcsulc_values',subject{1});
        %save(sprintf('%s_rh_calcsulc_values_v2.mat',subject{1}),'output');
        %calcSulc_save
   %end
%end

