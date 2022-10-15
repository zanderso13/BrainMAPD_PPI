function run_subject_firstlevel_BrainMAPD_rest(PID)
%% var set up
if nargin==0 % defaults just for testing 
    % Define some 
    PID = "10008";
    
end

are_you_doing_activation_first_levels = 1;
are_you_doing_ppi_first_levels = 0;

contrast = 'anticipation';
seed_region = 'OldhamAntRewVS'; % OldhamAntRewVS OldhamAntLossVS
overwrite = 1;
ses = 2;
run = 1;

% Define some paths
basedir = '/projects/b1108/projects/brainmapd_ppi';

% directories
% first is where your stats files will be output to
fl_dir = fullfile(basedir,'/first_levels/activation');
% next is where the preprocessed data is
preproc_dir = '/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid/smoothed_MID';
% where the raw data lives (raw meaning before preprocessing)
raw_dir = '/projects/b1108/data/BrainMAPD';
% where framewise displacement files will be saved
save_dir = fullfile(basedir,'/first_levels/FD');
% directory where I'm storing timing files for the MID
timing_dir = fullfile(strcat('/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid_spm_timing/run-',num2str(run)),contrast);


numPID = num2str(PID);
PID = strcat('sub-',numPID);


fprintf(['Preparing 1st level model for MID task for ' PID ' / ' ses], ['Overwrite = ' num2str(overwrite)]);


ndummies = 2;
TR = 2.05;

%% Model for MID task

% FL directory for saving 1st level results: beta images, SPM.mat, etc.
in{1} = {fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), 'mid', strcat('run-', num2str(run)))};

% preproc images
rundir = fullfile(preproc_dir, PID, strcat('ses-',num2str(ses)), 'func');
in{2} = cellstr(spm_select('ExtFPList', preproc_dir, strcat('^ssub-',numPID,'.*task-MID_run-',num2str(run),'_space-MNI152NLin6Asym_desc-preproc_bold.nii'), ndummies+1:9999));
    
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
if are_you_doing_activation_first_levels == 1
    % fmriprep output
    confound_fname = filenames(fullfile(preproc_dir, strcat(PID,'*task-MID*confounds*.tsv')));

    % find the raw image file, for spike detection
    raw_img_fname = filenames(fullfile(raw_dir, PID, strcat('ses-',num2str(ses)), 'func', strcat('*task-MID_run-0',num2str(run),'_bold.nii*')));
    cd(fullfile(preproc_dir));
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

    confoundFile = fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), 'mid', strcat('run-', num2str(run)), 'mid_confounds.mat');

    in{4} = {confoundFile};

    % checks
    if any(cellfun( @(x) isempty(x{1}), in))
        in
        error('Some input to the model is missing')
    end

end

if are_you_doing_ppi_first_levels == 1

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
    save(fullfile(fl_dir, PID, strcat('ses-',num2str(ses)), 'mid', strcat('run-', num2str(run)),'mid_confounds.mat'),'R','names');

    
    % run spm FL estimation
    cwd = pwd;
    job = strcat('MID_SPM_',contrast,'_template.m');
    %%
    spm('defaults', 'FMRI')
    spm_jobman('serial',job,'',in{:});

    cd(cwd);
end

end


