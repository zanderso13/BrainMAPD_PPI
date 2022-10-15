hyperalignment = 0;

load('/Volumes/DataCave/ACNlab/BrainMAPD/RS/outcomes/trilevel.mat')
load('/Volumes/DataCave/ACNlab/BrainMAPD/RS/outcomes/immune_data.mat')
load('/Volumes/DataCave/ACNlab/BrainMAPD/RS/outcomes/meds.mat')
load('/Volumes/DataCave/ACNlab/BrainMAPD/RS/outcomes/demographics.mat')
load('/Users/zaz3744/Documents/GitHub/BrainMAPD_resting_state/AAL3_1mm.mat')
load('/Volumes/DataCave/ACNlab/BrainMAPD/RS/outcomes/scan_day.mat')

exclusions = {'10001','10041','10034','10041',...
    '10059','10067','10081','10088','10111','10135',...
    '10141','10196','10282','10309','10336','10438','10439',...
    '20085','20108','20644','21238','21257','21268',...
    '20133','20460','20948','20309','20996','21001',...
    '21523','20309'};

fnames = filenames(fullfile('/Volumes/DataCave/ACNlab/BrainMAPD/RS/conn_matrices_ha/*hyp.mat'));

% remove mat files from subject list
for badsub = 1:length(exclusions)
    fnames(contains(fnames,exclusions{badsub}),:)=[];
end
    
nsubs = length(fnames);

for sub = 1:length(fnames)
    id=fnames{sub}(56:60);
    %currmotion = filenames(fullfile('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/FD',strcat('*',id,'*.mat')));
    %load(currmotion{1})
    % motion data
%     if isempty(FD)
%         allFD(sub,:)=NaN;
%     else
%         allFD(sub,:)=mean(FD.framewise_displacement);
%     end
    % symptom data
    if isempty(trilevel(trilevel.id==str2double(id),:))
        symp(sub,:)=NaN;
    else
        symp(sub,:)=trilevel(trilevel.id==str2double(id),:);
    end
    % immune data
    if isempty(immune.T1BDicsavg(immune.PID==str2double(id)))
        imm_comp(sub,1)=NaN;
    else
        imm_comp(sub,1)=immune.T1BDicsavg(immune.PID==str2double(id));
    end
    % med data
    if isempty(meds.T1SCcmipsychany(meds.PID==str2double(id)))
        curr_meds(sub,1)=0;
    else
        curr_meds(sub,1)=meds.T1SCcmipsychany(meds.PID==str2double(id));
    end
    % demographic data
    if str2double(id) < 20000
        site(sub,1) = 1;
    else
        site(sub,1) = 0;
    end
    % scan site
    if isempty(dem.sex(dem.PID==str2double(id)))
        sex(sub,1)=0;
        dob(sub,1)=NaN;
    else
        sex(sub,1)=dem.sex(dem.PID==str2double(id));
        dob(sub,1)=dem.dob(dem.PID==str2double(id));
    end

    % scan date
    if isempty(scan_day.T1M1date(scan_day.PID==str2double(id)))
        scanday(sub,1)=0;
    else
        scanday(sub,1)=scan_day.T1M1date(scan_day.PID==str2double(id));
    end
    temp = load(fnames{sub});
    all_data{:,:,sub} = temp.hyp_mat;
end

% load('300ROI_atlas_obj.mat')
% 
% SMD = strcmp(atl.labels,'SomatomotorDorsal');
% CO = strcmp(atl.labels,'CinguloOpercular');
% A = strcmp(atl.labels,'Auditory');
% DMN = strcmp(atl.labels,'DefaultMode');
% PM = strcmp(atl.labels,'ParietoMedial');
% V = strcmp(atl.labels,'Visual');
% FP = strcmp(atl.labels,'FrontoParietal');
% S = strcmp(atl.labels,'Salience');
% VA = strcmp(atl.labels,'VentralAttention');
% DA = strcmp(atl.labels,'DorsalAttention');
% MTL = strcmp(atl.labels,'MedialTemporalLobe');
% Rew = strcmp(atl.labels,'Reward');
% SML = strcmp(atl.labels,'SomatomotorLateral');

% input what network you want here. You can also do addition to look at
% several networks at once.
% seed to seed associations: 
% VS - mOFC: 
% idx1 = [246,247]; idx2 = [108,116];
% Amyg - mOFC: 
% idx1 = [244,245]; idx2 = [108,116];
% Amyg - VS: 
% idx1 = [246,247]; idx2 = [244,245];
% Network analysis example for connectivity within reward network: idx1 = Rew; idx2 = Rew;
% idx1 = Rew; idx2 = Rew;

% curr_analysis = all_data(idx1,idx2,:);
% 
% for sub = 1:size(curr_analysis,3)
%     conn(sub,1) = mean2(curr_analysis(:,:,sub));
% end
% 
% dat = [symp,array2table(imm_comp),array2table(conn),array2table(sex),array2table(curr_meds),array2table(site),array2table(allFD)]; 
% 
% fprintf(strcat('How many motion exclusions: ',num2str(sum(dat.allFD>0.3))))
% dat(dat.allFD>0.3,:)=[];
% 
% dat(isnan(dat.imm_comp),:) = [];
% if hyperalignment == 1
%     dat = dat(1:nsubs,:);
% end

%% age at first scan calculator

mean((scanday-dob)/8760); % 8760 is the number of hours in a year
std((scanday-dob)/8760);

%% hyperalignment
if hyperalignment == 1
    %PID_list = dat.id;
    %[alignedvs1,alignedvs2,alignedamyg1,alignedamyg2,alignedofc1,alignedofc2,transformvs1, transformvs2, transformamyg1, transformamyg2, transformofc1, transformofc2] = BrainMAPD_hyperalign(PID_list);
    [aligned, transforms] = hyperalign(all_data{:});
end

for sub = 1:length(all_data)
    if hyperalignment == 1
        X(sub,:) = [aligned{sub}(:,1)',aligned{sub}(:,2)',aligned{sub}(:,3)',...
            aligned{sub}(:,4)',aligned{sub}(:,5)',aligned{sub}(:,6)',aligned{sub}(:,7)',...
            aligned{sub}(:,8)'];
    else
        X(sub,:) = [all_data{sub}(:,1)',all_data{sub}(:,2)',all_data{sub}(:,3)',...
            all_data{sub}(:,4)',all_data{sub}(:,5)',all_data{sub}(:,6)',all_data{sub}(:,7)',...
            all_data{sub}(:,8)'];
    end
    
end

mdl = fitrlinear(X,symp.GenDis,'KFold',5,...
    'Learner','leastsquares','Regularization','ridge');
%% Final cleaning
% dat(dat.Anhedon(:) > 3,:) = []; dat(dat.Anhedon(:) < -3,:) = [];
% dat(dat.GenDis(:) > 3,:) = []; dat(dat.GenDis(:) < -3,:) = [];
% dat(dat.Fears(:) > 3,:) = []; dat(dat.Fears(:) < -3,:) = [];
% dat(zscore(dat.conn)>3,:) = []; dat(zscore(dat.conn)<-3,:) = [];
% dat(dat.imm_comp>3,:) = []; dat(dat.imm_comp<-3,:) = [];


%% compare conventional analyses with same thing on aligned data

% mdl1a = fitlm(dat,'GenDis ~ conn + curr_meds + sex + site');
% mdl1b = fitlm(dat,'GenDis ~ conn + imm_comp + imm_comp*conn + curr_meds + sex + site');
% mdl2a = fitlm(dat,'Anhedon ~ conn + curr_meds + sex + site');
% mdl2b = fitlm(dat,'Anhedon ~ conn + imm_comp + imm_comp*conn + curr_meds + sex + site');
% mdl3a = fitlm(dat,'Fears ~ conn + curr_meds + sex + site');
% mdl3b = fitlm(dat,'Fears ~ conn + imm_comp + imm_comp*conn + curr_meds + sex + site');
% mdl4 = fitlm(dat,'conn ~ GenDis + Anhedon + Fears + curr_meds + sex + site');
% plotAdjustedResponse(mdl4,'Anhedon')

% if hyperalignment == 1
%     
%     for sub = 1:length(alignedofc2)
%         amyg(:,:,sub) = [zscore(alignedamyg1{sub}),zscore(alignedamyg2{sub})];
%         vs(:,:,sub) = [zscore(alignedvs1{sub}),zscore(alignedvs2{sub})];
%         ofc(:,:,sub) = [zscore(alignedofc1{sub}),zscore(alignedofc2{sub})];
%         ha_conn(:,:,sub) = corr(vs(:,:,sub),ofc(:,:,sub));
%     end
%     
%     dat = [dat,array2table(ha_conn)];
% 
% %     ha_mdl1 = fitlm(dat,'ha_conn ~ imm_comp + GenDis + imm_comp*GenDis + curr_meds + sex + site')
% %     ha_mdl2 = fitlm(dat,'ha_conn ~ imm_comp + Anhedon + imm_comp*Anhedon + curr_meds + sex + site')
% %     ha_mdl3 = fitlm(dat,'ha_conn ~ imm_comp + Fears + imm_comp*Fears + curr_meds + sex + site')
%     ha_mdl4 = fitlm(dat,'ha_conn ~ GenDis + Anhedon + Fears + curr_meds + sex + site')
%     
% end
