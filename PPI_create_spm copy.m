%% The true PPI script
% You have all the pieces, more or less. Now you need to 
% 1. Extract time series data from preproc'ed BOLD data
% 2. Input those time series into a new .mat file, with the rest of the
% nuisance regressors from previous SPM.mat file
% 3. Rerun first levels! Which means just resubmit them but with a new
% template file that doesn't specificy design but rather just inputs the
% entire thing as a multiple regressors option. Hopefully it lets you do
% this.

% Input the number of contrasts you defined
number_of_events = 11;
%%%%%%% USER DEFINED %%%%%%%%%%
% Define some paths
% This is going to generate a first level script to be submitted to the
% cluster with each run. Where do you want all these .sh scripts saved?
scriptdir = '/projects/b1108/projects/BrainMAPD_func_conn/first_levels/PPI/quest_submission';

% Where are all your scripts saved for first levels? i.e. where is the
% acnlab_repo folder? Also where is spm12... you need spm

repodir = '~/repo';

% directories
% first is where your stats files will be output to
outputdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/PPI/regressors';
% next is where the preprocessed data is
preprocdir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/fmriprep';
% need a masks directory where you'll put binary .nii files
roidir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/roi';
% where are the first levels? More importantly where are the SPM.mat files
fldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/func_conn/first_levels/first_level_output';


% What run of your task are you looking at?
run = 1;
% What session appears in your raw filenames when in BIDS format?
ses = 2;
% Last thing is janky but bear with me. How long are your participant ID's?
% i.e. 10234 would correspond with a 5 for this variable
ID_length = 5;

%%%%%%% END USER DEFINED %%%%%%%%%%

% This will give you PID, don't change this.
foldernames = char(filenames(fullfile(preprocdir,'sub*/')));
sublist = foldernames(:,size(foldernames,2)-(ID_length):size(foldernames,2)-1);

sublist = string(sublist);
% I like this better than the overwrite option. This way, the list of
% subjects to be run will change each time I run the script. The overwrite
% option will still be helpful here to turn things on or off.


% Run/submit first level script

%cd(scriptdir)
for sub = 1:length(sublist)

    spm_file = filenames(fullfile(fldir,strcat('sub-',sublist{sub}),strcat('ses-',num2str(ses)),strcat('run-',num2str(run)),'MID','SPM.mat'));
    run1 = fmri_data(filenames(fullfile(preprocdir,strcat('sub-',sublist{sub}),strcat('ses-',num2str(ses)),'func','*MID*run-01*_bold.nii')));
    roi = fmri_data(filenames(fullfile(roidir,'*OFC*nii')));
    load(spm_file{1});

    time_series = run1.dat;
    roi_dat = roi.dat;
    roi_time_series = time_series .* roi_dat;

    mean_roi_bold = sum(roi_time_series,1) ./ 281;

    std_mean_roi_bold = zscore(mean_roi_bold');
    
    
    roi_time_regressor = std_mean_roi_bold(3:281);
    regressors = SPM.xX.X;


    % convolve task*HRF with Bold time series

    
    con_names = {'Loss Anticipation','Loss Anticipation Zero','Gain Anticipation','Gain Anticipation Zero','Loss Neutral','Loss Consumption','Gain Consumption','Gain Neutral','Motor'};
    con_names_for_save = {'Loss_Anticipation','Loss_Anticipation_Zero','Gain_Anticipation','Gain_Anticipation_Zero','Loss_Neutral','Loss_Consumption','Gain_Consumption','Gain_Neutral','Motor'};
    for con = 1:number_of_events
        con_temp(:,con) = double(regressors(:,con) .* roi_time_regressor);        
    end
    
    figure();
    plot(con_temp)
    title('Contrasts * Preprocessed BOLD time course data')
    temp_plot_name = strcat(sublist{sub},'_','contrasts.jpg');
        
    % save the figure if you'd like
    %saveas(gcf,fullfile(save_dir,temp_plot_name))
    % save the regressor file
    %close all
    % concatenate all the regressors
    R = [con_temp,regressors(:,12:size(regressors,2))];
    

    
    temp_file_name = fullfile(outputdir,strcat(sublist{sub},'_PPI.mat'));
    save(temp_file_name,'R')
end

