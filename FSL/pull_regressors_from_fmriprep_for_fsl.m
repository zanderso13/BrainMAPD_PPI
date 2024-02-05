% note: matlab doesn't like tsv's, so I've renamed all files to have .txt
% extension. This makes it so that the readtable function works
% appropriately.

confounddir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/globus_transfer_confound_files';
outputdir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/first_level_confounds_dissertation';

cd(confounddir)
id_length = 9; % length of sub id. Ex: sub-10001
% drop 10 images for rest and 2 for mid
number_of_images_dropped = 2; % fsl wants the extra confounds file to match time course data after it removes images

fnames = filenames(fullfile('sub-*/ses-2/func/*MID_run-2*.txt'));


for sub = 1:length(fnames)
    id{sub,1} = fnames{sub}(1:id_length);
    T = readtable(fnames{sub});

    GSR = T.global_signal;
    FD(sub,1) = nanmean(T.framewise_displacement);
    DVARS = T.dvars;
    
    CSF = T(:,contains(T.Properties.VariableNames,'csf'));
    WM = T(:,contains(T.Properties.VariableNames,'white_matter'));
    outliers = T(:,contains(T.Properties.VariableNames,'motion'));
    transx = T(:,contains(T.Properties.VariableNames,'trans_x'));
    transy = T(:,contains(T.Properties.VariableNames,'trans_y'));
    transz = T(:,contains(T.Properties.VariableNames,'trans_z'));
    rotx = T(:,contains(T.Properties.VariableNames,'rot_x'));
    roty = T(:,contains(T.Properties.VariableNames,'rot_y'));
    rotz = T(:,contains(T.Properties.VariableNames,'rot_z'));
    
    % 6 motion regressors
    % totmotion = [transx(:,1),transy(:,1),transz(:,1),rotx(:,1),roty(:,1),rotz(:,1),outliers]; %totmotion = table2array(totmotion);
    % 24 motion regressors
    totmotion = [transx,transy,transz,rotx,roty,rotz,outliers]; %totmotion = table2array(totmotion);


    R = [totmotion,WM(:,1:4),CSF(:,1:4)];
    names = R.Properties.VariableNames;
    R = table2array(R);
    
    R(isnan(R)) = 0;
    R = R(number_of_images_dropped+1:size(T,1),:);
    confounds_fname = strcat(id{sub},'_rest_ses-2_run-1_confounds.mat'); % _mid_ses-2_run-1_confounds.mat
%     motionfname = strcat(id,'motion.txt');
%     GSRfname = strcat(id,'_gsr.txt');
%     WMfname = strcat(id,'_wm.txt');
%     CSFfname = strcat(id,'_csf.txt');
%     FDfname = strcat(id,'_fd.txt');
%     DVARSfname = strcat(id,'_dvars.txt');
    
%     if ~exist(fullfile(outputdir,id{sub}))
%         mkdir(fullfile(outputdir,id{sub}))
%     end
%     savelocation = fullfile(outputdir,confounds_fname);
%     save(savelocation, 'R','names')
%     writematrix(GSR,fullfile(outputdir,id,GSRfname))
%     writematrix(WM,fullfile(outputdir,id,WMfname))
%     writematrix(CSF,fullfile(outputdir,id,CSFfname))
%     writematrix(FD,fullfile(outputdir,id,FDfname))
%     writematrix(DVARS,fullfile(outputdir,id,DVARSfname))

end

id = cell2table(id); FD = array2table(FD); motion_table = [id,FD];
