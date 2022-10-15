% run first level models for REST task for a subject
% the URSI and whether to overwrite existing first levels.  Overwrite = 0
% (default) or 1. before running this file, must run
% read_timings_make_onsets.m

    
 
function run_subject_firstlevel_BrainMAPD_rest(PID)
%% var set up
if nargin==0 % defaults just for testing 
    % Define some 
    PID = "10006";
    
end

contrast = 'anticipation';
seed_region = 'OldhamAntRewVS'; % OldhamAntRewVS OldhamAntLossVS
overwrite = 1;
ses = 2;
run = 1;

% Define some paths
basedir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';

% directories
% first is where your stats files will be output to
directories{1} = fullfile(basedir,'/first_levels/activation');
% next is where the preprocessed data is
directories{2} = fullfile(basedir,'fmriprep');
% where the raw data lives (raw meaning before preprocessing)
directories{3} = '/Volumes/DataCave/ACNlab/BrainMAPD/MID_2020/projects/b1108/data/BrainMAPD';
% directories{3} = '/home/zaz3744/ACNlab/data/MWMH';
% where framewise displacement files will be saved
directories{4} = fullfile(basedir,'/first_levels/FD');


fl_dir = directories{1};
preproc_dir = directories{2};
raw_dir = directories{3};
save_dir = directories{4};

numPID = num2str(PID);
PID = strcat('sub-',numPID);

timing_dir = fullfile(strcat('/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/mid_spm_timing/run-',num2str(run)),contrast);

fprintf(['Preparing 1st level model for MID task for ' PID ' / ' ses], ['Overwrite = ' num2str(overwrite)]);


ndummies = 2;
TR = 2.05;

%% Model for REST task

% FL directory for saving 1st level results: beta images, SPM.mat, etc.
in{1} = {fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), strcat('run-', num2str(run)), 'mid')};

% preproc images
rundir = fullfile(preproc_dir, PID, strcat('ses-',num2str(ses)), 'func');
in{2} = cellstr(spm_select('ExtFPList', rundir, strcat('^ssub.*task-MID_run-',num2str(run),'_space-MNI152NLin6Asym_desc-preproc_bold.nii'), ndummies+1:9999));
    
if isempty(in{2}{1})
    warning('No preprocd functional found')
    return
end

% onset files
in{3} = filenames(fullfile(timing_dir, strcat(numPID,'*')));

if isempty(in{3})
    warning('No modeling found (behav data might be missing)')
    return
end

%% nuisance covs

% fmriprep output
confound_fname = filenames(fullfile(preproc_dir, strcat(PID), strcat('ses-',num2str(ses)), 'func', '*task-MID*confounds*.tsv'));

% find the raw image file, for spike detection
raw_img_fname = filenames(fullfile(raw_dir, PID, strcat('ses-',num2str(ses)), 'func', strcat('*task-MID_run-0',num2str(run),'_bold.nii*')));
cd(fullfile(preproc_dir, PID));
confound_savedir = fullfile(fl_dir, PID, strcat('confounds-MID_run-0',num2str(run),'.txt'));
% get nuis covs

[Rfull, Rselected, n_spike_regs, framewise_displacement_final, gsr_final] = make_nuisance_for_fl(confound_fname{run}, raw_img_fname, TR, 0, confound_savedir, preproc_dir);
save(fullfile(save_dir, strcat(PID, '_ses', num2str(ses), '_run', num2str(run), '.mat')), 'framewise_displacement_final')


% choose which matrix to use
R = Rselected;

% its now possible that some of the spike regs are all zero, b/c the spikes
% were discarded in the step above. find all-zero regs and drop
R(:, ~any(table2array(R))) = [];
R = R(ndummies+1:end, :); %discard dummy vols

% put in SPM format: matrix called 'R', and 'names'
names = R.Properties.VariableNames;
R = table2array(R);

confoundFile = fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), strcat('run-', num2str(run)),'mid','mid_confounds.mat');

in{4} = {confoundFile};

% checks
if any(cellfun( @(x) isempty(x{1}), in))
    in
    error('Some input to the model is missing')
end


% check for SPM.mat and overwrite if needed
skip = 0;
if exist(fullfile(in{1}{1},'SPM.mat'),'file')
    if overwrite
        fprintf('\n\nWARNING: EXISTING SPM.MATAND BETA FILES WILL BE OVERWRITTEN\n%s\n\n',fullfile(in{1}{1},'SPM.mat'));
        rmdir(in{1}{1},'s');
    else
        fprintf('\n\nFirst levels already exist, wont ovewrite: %s\n\n',fullfile(in{1}{1},'SPM.mat'));
        skip=1;
    end
end

if ~skip
    % make dir for beta and contrast files
    if ~isdir(in{1}{1}), mkdir(in{1}{1}); end
    save(fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), strcat('run-', num2str(run)), 'mid','mid_confounds.mat'),'R','names');


    % run spm FL estimation
    cwd = pwd;
    job = 'MID_SPM_anticipation_template.m';
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job,'',in{:});

    cd(cwd);
end

end


