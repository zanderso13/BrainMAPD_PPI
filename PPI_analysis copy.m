trial_type = 'con'; % possible input: 'ant' 'con'
specificity = 1; % always keep this on. No specificity is outputting something a little strange right now
perm_test = 0;
plots = 1;
% Exclusions:I'm removing overlap manually

%corrupt_data = {'10081','10327','20309','20235'};
% 0.3 cutoff exclusion
% motion_problems = {'10074','10111','10161','10196','10236','10264','10272','10274','10282','10296','10308','10341','10434','10481','21597','21684'};
% 
% admin_problems = {'10002','10029','10141','10230','10272','10307','10332','20025','20085','20123','20143','20150','20302','20456','20610','20652','20655','20877','20936','21020','21081','21126'};
% 
% grubbs_outlier = {'21386'}; % outlier based on anhedonia scores


basedir = '/Volumes/DataCave/ACNlab/BrainMAPD/PPI/PPI/ppi_fldir/consumption';
maskdir = '/Users/zaz3744/Documents/current_projects/ACNlab/masks/ROI_BrainMAPD_functional/consumption';
% roidir = {'OldhamOFC','OldhamAntRewVS','HOAmyg'};
roidir = {'OldhamConsumpVS'}; % 'OldhamAntRewVS' 'OldhamConsumpVS'

%exclusions = [corrupt_data,motion_problems,admin_problems,grubbs_outlier];
 exclusions = {'10001','10084','10094','10102','10125','10140'...
    '10143','10148','10161','10264','10274','10296','10319',...
    '10327','10422','10423','10434','10443','10461','10471',...
    '20032','20050','20108','20309','20317','20507','20564',...
    '20674','21111','21178','21223','21597','21386'};
% load masks in this section
if strcmp(trial_type,'ant')
    mOFC = fmri_data(filenames(fullfile(maskdir,'OFC*Oldham*nii')));
    VSloss = fmri_data(filenames(fullfile(maskdir,'VS*Oldham_Loss*nii')));
    VSwin = fmri_data(filenames(fullfile(maskdir,'VS*Oldham_Rew*nii')));
    amyg = fmri_data(filenames(fullfile(maskdir,'HO_Amygdala*nii')));
else
    mOFC = fmri_data(filenames(fullfile(maskdir,'OFC*Oldham*nii')));
    amyg = fmri_data(filenames(fullfile(maskdir,'HO_Amygdala*nii')));
    VSconsump = fmri_data(filenames(fullfile(maskdir,'VS*Oldham_Con.nii')));
end
% load regressors
clinicaldir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/clinical_data';
med_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';
demo_dir = '/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD';
load(fullfile(clinicaldir,'first_year_project_trilevel_T1.mat'));
load(fullfile(med_dir,'Medication_T2.mat'));
load(fullfile(demo_dir,'demographics.mat'));
load(fullfile(clinicaldir,'first_year_project_BrainMAPD_clinical_diagnoses_final.mat'))
load(fullfile(demo_dir,'SES.mat'))
for roi = 1:length(roidir)
    
    % Pull out brain data for all included subjects
    cd(fullfile(basedir,roidir{roi}))
    
    con1run1 = filenames(fullfile('*/ses-2/run-1/MID/con_0001.nii'));
    con2run1 = filenames(fullfile('*/ses-2/run-1/MID/con_0002.nii'));
    
    con1run2 = filenames(fullfile('*/ses-2/run-2/MID/con_0001.nii'));
    con2run2 = filenames(fullfile('*/ses-2/run-2/MID/con_0002.nii'));
    
    % for run1
    for ex = 1:length(exclusions)
        ind1 = contains(con1run1,exclusions(ex));
        con1run1(ind1) = [];
        ind2 = contains(con2run1,exclusions(ex));
        con2run1(ind2) = [];
    end
    % for run2
    for ex = 1:length(exclusions)
        ind1 = contains(con1run2,exclusions(ex));
        con1run2(ind1) = [];
        ind2 = contains(con2run2,exclusions(ex));
        con2run2(ind2) = [];
    end
    % match across runs
    for sub = 1:length(con1run2)
        curr_sub = con1run2{sub}(1:5);
        curr_ind = contains(con1run1,curr_sub);
        if sum(curr_ind) == 1
            tempbeta1{sub} = con1run1{curr_ind};
            tempbeta2{sub} = con2run1{curr_ind};
        else
            disp(con1run2{sub}(1:5))
            continue
        end
    end

    con1run1 = tempbeta1';
    con2run1 = tempbeta2';
    
    % load in data after removing exclusions
    % run 1
    con1run1 = fmri_data(con1run1);  
    con2run1 = fmri_data(con2run1);
   
    % run 2
    con1run2 = fmri_data(con1run2);
    con2run2 = fmri_data(con2run2);
    
    % average the runs
    acon1.(roidir{roi}) = con1run1; acon1.(roidir{roi}).dat = acon1.(roidir{roi}).dat+con1run2.dat; acon1.(roidir{roi}).dat = acon1.(roidir{roi}).dat .* 0.5;
    acon2.(roidir{roi}) = con2run1; acon2.(roidir{roi}).dat = acon2.(roidir{roi}).dat+con2run2.dat; acon2.(roidir{roi}).dat = acon2.(roidir{roi}).dat .* 0.5;
    
    % extract region of interest estimates for seed to seed analysis
    if trial_type == 'ant'
        % con 1
        c1r1amyg = extract_roi_averages(acon1.(roidir{roi}),amyg);
        c1r1VSloss = extract_roi_averages(acon1.(roidir{roi}),VSloss);
        c1r1VSrew = extract_roi_averages(acon1.(roidir{roi}),VSwin);
        c1r1mOFC = extract_roi_averages(acon1.(roidir{roi}),mOFC);
        % con 2
        c2r1amyg = extract_roi_averages(acon2.(roidir{roi}),amyg);
        c2r1VSloss = extract_roi_averages(acon2.(roidir{roi}),VSloss);
        c2r1VSrew = extract_roi_averages(acon2.(roidir{roi}),VSwin);
        c2r1mOFC = extract_roi_averages(acon2.(roidir{roi}),mOFC);
    else
        % con 1
        c1r1amyg = extract_roi_averages(acon1.(roidir{roi}),amyg);
        c1r1VSconsump = extract_roi_averages(acon1.(roidir{roi}),VSconsump);
        c1r1mOFC = extract_roi_averages(acon1.(roidir{roi}),mOFC);
        % con 2
        c2r1amyg = extract_roi_averages(acon2.(roidir{roi}),amyg);
        c2r1VSconsump = extract_roi_averages(acon2.(roidir{roi}),VSconsump);
        c2r1mOFC = extract_roi_averages(acon2.(roidir{roi}),mOFC);
    end
    
    
    % Extract clinical, demographic, medication info
    
    trilevel_array = [trilevel_T1.id,trilevel_T1.GenDis,trilevel_T1.Anhedon,trilevel_T1.Fears];
    med_array = [T2MedicationInventory3(:,1),T2MedicationInventory3(:,88)];
    med_array = table2array(med_array);
    dem_array = [BrainMAPDT1S1Demo.PID,BrainMAPDT1S1Demo.sex,BrainMAPDT1S1Demo.ethnicity,BrainMAPDT1S1Demo.race,BrainMAPDT1S1Demo.race2];
    dsm_array = table2array(clinical_info);
    
    % pull data for real subjects
    for sub = 1:length(tempbeta1)
        % trilevel 
        PID = str2num(tempbeta1{sub}(1:5));
        if isempty(find(trilevel_T1.id(:) == PID)) == 0
            curr = find(trilevel_T1.id(:) == PID);
            trilevel_regressors(sub,:) = trilevel_array(curr,:);
        else
            % for OFC loss_ppi_fnames{sub}(100:104); for HO_VMPFC loss_ppi_fnames{sub}(105:109)
            disp(strcat(num2str(PID), ' missing clinical info')) 
            trilevel_regressors(sub,:) = NaN;
            trilevel_regressors(sub,1) = PID;
        end
        % dsm
        if isempty(find(clinical_info.PID(:) == PID)) == 0
            curr = find(clinical_info.PID(:) == PID);
            dsm_regressors(sub,:) = dsm_array(curr,2:4);
        else
            % for OFC loss_ppi_fnames{sub}(100:104); for HO_VMPFC loss_ppi_fnames{sub}(105:109)
            disp(strcat(num2str(PID), ' missing clinical info')) 
            dsm_regressors(sub,:) = [0,0,0];
            %dsm_regressors(sub,1) = PID;
        end
        % medication
        if isempty(find(med_array(:,1) == PID)) == 0
            curr2 = find(med_array(:,1) == PID);
            med_regressors(sub,:) = med_array(curr2,2);
        else
            disp(strcat(num2str(PID), ' missing medication info'))
            med_regressors(sub,:) = [0];
        end
        % demographics
        if isempty(find(dem_array(:,1) == PID)) == 0
            curr2 = find(dem_array(:,1) == PID);
            demographic_regressors(sub,:) = dem_array(curr2,2);%:5);
        else
            disp(strcat(num2str(PID), ' missing demographic info'))
            demographic_regressors(sub,:) = 0;
        end
        % site 
        if PID < 20000
            site_regressor(sub,1) = 1;
        else
            site_regressor(sub,1) = 0;
        end
        % SES
        if isempty(find(SES.PID == PID)) == 0
            curr_job(sub) = SES.T1SESwsstatus(find(SES.PID == PID));
            curr_inc(sub) = SES.T1SESfaminc(find(SES.PID == PID));
        else
            curr_job(sub) = 0;
            curr_inc(sub) = NaN;
        end
    end
    % pull demographic and other data for excluded subs
    for sub = 1:length(exclusions)
        % trilevel 
        PID = str2num(exclusions{sub}(1:5));
        if isempty(find(trilevel_T1.id(:) == PID)) == 0
            curr = find(trilevel_T1.id(:) == PID);
            ex_trilevel_regressors(sub,:) = trilevel_array(curr,:);
        else
            % for OFC loss_ppi_fnames{sub}(100:104); for HO_VMPFC loss_ppi_fnames{sub}(105:109)
            disp(strcat(num2str(PID), ' missing clinical info')) 
            ex_trilevel_regressors(sub,:) = NaN;
            ex_trilevel_regressors(sub,1) = PID;
        end
        % dsm
        if isempty(find(clinical_info.PID(:) == PID)) == 0
            curr = find(clinical_info.PID(:) == PID);
            ex_dsm_regressors(sub,:) = dsm_array(curr,2:4);
        else
            % for OFC loss_ppi_fnames{sub}(100:104); for HO_VMPFC loss_ppi_fnames{sub}(105:109)
            disp(strcat(num2str(PID), ' missing clinical info')) 
            ex_dsm_regressors(sub,:) = [0,0,0];
            %dsm_regressors(sub,1) = PID;
        end
        % medication
        if isempty(find(med_array(:,1) == PID)) == 0
            curr2 = find(med_array(:,1) == PID);
            ex_med_regressors(sub,:) = med_array(curr2,2);
        else
            disp(strcat(num2str(PID), ' missing medication info'))
            ex_med_regressors(sub,:) = [0];
        end
        % demographics
        if isempty(find(dem_array(:,1) == PID)) == 0
            curr2 = find(dem_array(:,1) == PID);
            ex_demographic_regressors(sub,:) = dem_array(curr2,2);%:5);
        else
            disp(strcat(num2str(PID), ' missing demographic info'))
            ex_demographic_regressors(sub,:) = 0;
        end
        % site 
        if PID < 20000
            ex_site_regressor(sub,1) = 1;
        else
            ex_site_regressor(sub,1) = 0;
        end
    end
    
    % whole brain regression
    
%     trilevel_regressors = array2table(trilevel_regressors); trilevel_regressors.Properties.VariableNames = {'PID','gendis','anhed','fears'};
%     med_regressors = array2table(med_regressors); med_regressors.Properties.VariableNames = {'meds'};
%     demographic_regressors = array2table(demographic_regressors); demographic_regressors.Properties.VariableNames = {'sex'};
%     site_regressor = array2table(site_regressor); site_regressor.Properties.VariableNames = {'site'};
    if specificity == 1
        whole_brain_regressors = [trilevel_regressors(:,2:size(trilevel_regressors,2)),med_regressors,demographic_regressors,site_regressor];
    else
        whole_brain_regressors = [trilevel_regressors(:,symptom_to_analyze),med_regressors,demographic_regressors,site_regressor];
    end
    
    acon1.(roidir{roi}).X = whole_brain_regressors(:,2:size(whole_brain_regressors,2));
    acon2.(roidir{roi}).X = whole_brain_regressors(:,2:size(whole_brain_regressors,2));
    regout1.(roidir{roi}) = regress(acon1.(roidir{roi}),[0.05, 'fdr']);
    regout2.(roidir{roi}) = regress(acon2.(roidir{roi}),[0.05, 'fdr']);
    
    % region of interest analysis
    roi_regressors_ofc = [whole_brain_regressors,c2r1mOFC.dat];
    roi_regressors_amyg = [whole_brain_regressors,c2r1amyg.dat];
    if trial_type == 'ant'
        roi_regressors_vs = [whole_brain_regressors,c2r1VSrew.dat];
    else
        roi_regressors_vs = [whole_brain_regressors,c2r1VSconsump.dat];
    end

    trilevel_regressors = array2table(trilevel_regressors); trilevel_regressors.Properties.VariableNames = {'PID','gendis','anhed','fears'};
    med_regressors = array2table(med_regressors); med_regressors.Properties.VariableNames = {'meds'};
    demographic_regressors = array2table(demographic_regressors); demographic_regressors.Properties.VariableNames = {'sex'};
    site_regressor = array2table(site_regressor); site_regressor.Properties.VariableNames = {'site'};
    amyg = array2table(c1r1amyg.dat); amyg.Properties.VariableNames = {'predict_amyg'};
    amyg_neg = array2table(c2r1amyg.dat); amyg_neg.Properties.VariableNames = {'control_amyg'};
    ofc = array2table(c1r1mOFC.dat); ofc.Properties.VariableNames = {'predict_ofc'};
    ofc_neg = array2table(c2r1mOFC.dat); ofc_neg.Properties.VariableNames = {'control_ofc'};
    s2s = [trilevel_regressors,med_regressors,demographic_regressors,site_regressor,amyg,ofc,amyg_neg,ofc_neg];
    
    
    if strcmp(roidir{roi}, 'OldhamAntRewVS') == 1
        mdl1 = fitlm(s2s, 'predict_ofc ~ gendis + anhed + fears + site + sex + meds + control_ofc')
        mdl2 = fitlm(s2s, 'predict_amyg ~ gendis + anhed + fears + site + sex + meds + control_amyg')
        mdl3 = fitlm(s2s, 'predict_ofc ~ gendis + site + sex + meds + control_ofc')
        mdl4 = fitlm(s2s, 'predict_ofc ~ anhed + site + sex + meds + control_ofc')
        mdl5 = fitlm(s2s, 'predict_ofc ~ fears + site + sex + meds + control_ofc')
        mdl6 = fitlm(s2s, 'predict_amyg ~ gendis + site + sex + meds + control_amyg')
        mdl7 = fitlm(s2s, 'predict_amyg ~ anhed + site + sex + meds + control_amyg')
        mdl8 = fitlm(s2s, 'predict_amyg ~ fears + site + sex + meds + control_amyg')
    elseif strcmp(roidir{roi}, 'OldhamConsumpVS') == 1
        mdl1 = fitlm(s2s, 'predict_ofc ~ gendis + anhed + fears + site + sex + meds + control_ofc')
        mdl2 = fitlm(s2s, 'predict_amyg ~ gendis + anhed + fears + site + sex + meds + control_amyg')
        mdl3 = fitlm(s2s, 'predict_ofc ~ gendis + site + sex + meds + control_ofc')
        mdl4 = fitlm(s2s, 'predict_ofc ~ anhed + site + sex + meds + control_ofc')
        mdl5 = fitlm(s2s, 'predict_ofc ~ fears + site + sex + meds + control_ofc')
        mdl6 = fitlm(s2s, 'predict_amyg ~ gendis + site + sex + meds + control_amyg')
        mdl7 = fitlm(s2s, 'predict_amyg ~ anhed + site + sex + meds + control_amyg')
        mdl8 = fitlm(s2s, 'predict_amyg ~ fears + site + sex + meds + control_amyg')
    end
    
    % Name variables
%     if specificity == 1
%         mdl1.PredictorNames{1} = 'General Distress';
%         mdl1.PredictorNames{2} = 'Anhedonia';
%         mdl1.PredictorNames{3} = 'Fears';
%         mdl1.PredictorNames{4} = 'Medication Status';
%         mdl1.PredictorNames{5} = 'Sex';
%         mdl1.PredictorNames{6} = 'Site';
%         mdl1.PredictorNames{7} = 'Loss Condition';
%         
%         mdl2.PredictorNames{1} = 'General Distress';
%         mdl2.PredictorNames{2} = 'Anhedonia';
%         mdl2.PredictorNames{3} = 'Fears';
%         mdl2.PredictorNames{4} = 'Medication Status';
%         mdl2.PredictorNames{5} = 'Sex';
%         mdl2.PredictorNames{6} = 'Site';
%         mdl2.PredictorNames{7} = 'Loss Condition';
%     end
end
    

%% permutation testing 
% This part of the script isn't dynamic. I've already determined the
% significant results and now I'm hard coding in the perm testing. Would be
% cool to write a different version of this someday that reads in
% significance stats and automatically initiates perm testing if needed

if perm_test == 1
    s2s_temp = s2s;
    for iter = 1:10000
        % create random indices
        idx = randperm(length(s2s_temp.predict_amyg));
        s2s_temp.predict_amyg = s2s.predict_amyg(idx);
        s2s_temp.predict_ofc = s2s.predict_ofc(idx);
        % generate temporary models
        tmdl3 = fitlm(s2s_temp, 'predict_ofc ~ gendis + site + sex + meds + control_ofc');
        tmdl4 = fitlm(s2s_temp, 'predict_ofc ~ anhed + site + sex + meds + control_ofc');
        tmdl5 = fitlm(s2s_temp, 'predict_ofc ~ fears + site + sex + meds + control_ofc');
        tmdl6 = fitlm(s2s_temp, 'predict_amyg ~ gendis + site + sex + meds + control_amyg');
        tmdl7 = fitlm(s2s_temp, 'predict_amyg ~ anhed + site + sex + meds + control_amyg');
        tmdl8 = fitlm(s2s_temp, 'predict_amyg ~ fears + site + sex + meds + control_amyg');
        
        % This is gonna be a lot of stuff, so I'll only grab the models
        % that are significant. If you ever have to do this for additional
        % measures, just copy and paste

        % VS-mOFC
        null_dist_ofc_gendis(iter) = tmdl3.Coefficients.Estimate(2);
        null_dist_ofc_anhed(iter) = tmdl4.Coefficients.Estimate(2);
        % VS-amyg
        null_dist_amyg_gendis(iter) = tmdl6.Coefficients.Estimate(2);
        null_dist_amyg_anhed(iter) = tmdl7.Coefficients.Estimate(2);
            
    end
    
    q_ofc_gendis = sum(null_dist_ofc_gendis < mdl3.Coefficients.Estimate(2)) / 10000;
    q_ofc_anhed = sum(null_dist_ofc_anhed > mdl4.Coefficients.Estimate(2)) / 10000;
    q_amyg_anhed = sum(null_dist_amyg_anhed > mdl7.Coefficients.Estimate(2)) / 10000;
    q_amyg_gendis = sum(null_dist_amyg_gendis < mdl6.Coefficients.Estimate(2)) / 10000;

    
end
    
%% plotting for publication
if plots == 1
    if specificity == 1
        fig1 = figure(); 
        title('General Distress')
        subplot(1,2,1); % vs-mofc
        plotAdjustedResponse(mdl1,'gendis','Marker','.','LineWidth',1); xlabel(' '); ylabel(' '); title(' '); legend off
        subplot(1,2,2); % vs-amygdala
        plotAdjustedResponse(mdl2,'gendis','Marker','.','LineWidth',1); xlabel(' '); ylabel(' '); title(' '); legend off
        
        fig1 = figure(); 
        title('Anhedonia-Apprehension')
        subplot(1,2,1); % vs-mofc
        plotAdjustedResponse(mdl1,'anhed','Marker','.','LineWidth',1); xlabel(' '); ylabel(' '); title(' '); legend off
        subplot(1,2,2); % vs-amygdala
        plotAdjustedResponse(mdl2,'anhed','Marker','.','LineWidth',1); xlabel(' '); ylabel(' '); title(' '); legend off
        
        fig1 = figure(); 
        title('Fears')
        subplot(1,2,1); % vs-mofc
        plotAdjustedResponse(mdl1,'fears','Marker','.','LineWidth',1); xlabel(' '); ylabel(' '); title(' '); legend off
        subplot(1,2,2); % vs-amygdala
        plotAdjustedResponse(mdl2,'fears','Marker','.','LineWidth',1); xlabel(' '); ylabel(' '); title(' '); legend off
        
    end
end
    
        
        
        
        