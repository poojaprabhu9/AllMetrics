% GLM fit including age  for each metrics

clear 
close all
%%
%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('C:/Users/pprabhu/OneDrive - UCSF/tresting/toolbox'))
addpath(genpath('C:/Users/pprabhu/OneDrive - UCSF/tresting/mfiles_forPooja/'))
n_reg=40;

common_source_path='C:/Users/pprabhu/OneDrive - UCSF/tresting/research_meg12/kranasinghe/mfiles/mfiles_from_Kudo5_60s/meganalysis/Database/Standard_mri_grid_vol';
atlasinfo = load(fullfile(common_source_path, 'atlasinfo'));
%%
outpath='C:/Users/pprabhu/OneDrive - UCSF/tresting/datafiles_PSD_dPTE/outputs_dell/Region_level_stats/sig_values/';
figpath='C:/Users/pprabhu/OneDrive - UCSF/tresting/datafiles_PSD_dPTE/outputs_dell/Region_level_stats/Figures';
datapath='C:/Users/pprabhu/OneDrive - UCSF/tresting/datafiles_PSD_dPTE/';
%%
p_thresh=0.05;
metrics='relative_pow';
freq={'deltatheta','alpha','beta'};
%% Plotting parameters
plttype = 'full';
adjustCM = 1;
atlasname='BNmod';
% these values are set intuitively looking at the results and by Kamalini's
% ppt
%% picking datafiles
AD=load(fullfile(datapath,'all_metrics_AD'));
CONT=load(fullfile(datapath,'all_metrics_CONT'));

% these values are set intuitively looking at the results and by Kamalini's
% ppt
rangeVal.deltatheta.AD = [15 40] ;
rangeVal.deltatheta.ADCONT= [-22 22];
rangeVal.alpha.AD = [10 30] ;
rangeVal.alpha.ADCONT= [-10 10];
rangeVal.beta.AD = [20 45] ;
rangeVal.beta.ADCONT= [-10 10];

%% dataSheets for demographic details
ADsheet = readtable(fullfile(datapath,"RadID_77AD_asPer_processed_data.xlsx"));
CONTsheet = readtable(fullfile(datapath,"RadID_90CONT_asPer_processed_data.xlsx"));

AD_age = ADsheet.Age_at_MEG;    
numericGender = cellfun(@(x) strcmp(x, 'Female'), ADsheet.Gender);  %female-1, male 0
AD_gender = double(numericGender);

CONT_age = CONTsheet.Age_at_MEG;    
numericGender = cellfun(@(x) strcmp(x, 'Female'), CONTsheet.Gender); 
CONT_gender = double(numericGender);

% group labels (1 for CONT and 2 for AD)
group_labels = [1 * ones(length(CONT_age), 1); 2 * ones(length(AD_age), 1)];

age = cat(1,CONT_age, AD_age);
gender = cat(1,CONT_gender, AD_gender);

nModules = [1:20 25:44];
%%
% results to store fitglm output
sig_vals = struct();
for f=1:length(freq)
    odir=fullfile(figpath,metrics,freq{f});
    %%%%%%%%%%% AD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fc_data_AD = AD.(metrics).modular.(freq{f});
    %%%%%%%%%%% Controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fc_data_CONT = CONT.(metrics).modular.(freq{f});
    %%%%%%%%%%%% final functional connectivity data %%%%%%%%%%%%%%%
    fc_data = cat(1,  fc_data_CONT, fc_data_AD);

    for reg = 1:n_reg
        % dependent (or response) variable
        y=fc_data(:,reg);

        % covariates + dependent (or response) variable in table
        T = table(age, gender, group_labels, y);
    
        % GLM analysis (assuming normal distribution)
        model_fit = fitglm(T, 'y ~ age + gender + group_labels'); % Fit linear model with all covariates
        
        % store p-values and t-values
        sig_vals.(metrics).(freq{f}).tval(1,reg) = model_fit.Coefficients.tStat('group_labels');
        sig_vals.(metrics).(freq{f}).pval(1,reg) = model_fit.Coefficients.pValue('group_labels');
        sig_vals.(metrics).(freq{f}).Coefficients{reg} = model_fit.Coefficients; 
       
    end
    % FDR correction
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(sig_vals.(metrics).(freq{f}).pval, p_thresh,'pdep','yes');
    corrected_pval_tmp=adj_p.*h;
    survived_tval_tmp=sig_vals.(metrics).(freq{f}).tval.*h;
    survived_tval = [survived_tval_tmp(1,1:20) 0 0 0 0 survived_tval_tmp(1,21:40) 0 0 0 0];
    corrected_pval = [corrected_pval_tmp(1,1:20) 0 0 0 0 corrected_pval_tmp(1,21:40) 0 0 0 0];
    
    afterPoint = extractAfter(num2str(p_thresh),'.'); 
    sig_vals.(metrics).(freq{f}).(['FDR' afterPoint 'tval']) = survived_tval;
    sig_vals.(metrics).(freq{f}).(['FDR' afterPoint 'pval']) = corrected_pval;

    %Visualisation %%%%%%
    % only 40 regions no subcortical regions
    avg_AD_data_tmp = squeeze(mean(fc_data_AD,1));
    avg_CONT_data_tmp = squeeze(mean(fc_data_CONT,1));

    % just making 48 regions for plotting
    avg_AD_data = [avg_AD_data_tmp(1,1:20) 0 0 0 0 avg_AD_data_tmp(1,21:40) 0 0 0 0];
    avg_CONT_data = [avg_CONT_data_tmp(1,1:20) 0 0 0 0 avg_CONT_data_tmp(1,21:40) 0 0 0 0];

     % plotting Control metric
    outputfilename=sprintf('avgMod_CONT_%s_%s',metrics,freq{f});
    range = [];
    range.minV = rangeVal.(freq{f}).AD(1); range.maxV = rangeVal.(freq{f}).AD(2);
    range.eps = 0; range.type = 2;
    [ctk]=[(othercolor('Mdarkrainbow',1000))]; CT = ctk;
    BNVplot(avg_CONT_data,range,plttype,CT,odir,atlasname,outputfilename,adjustCM) 

    % plotting AD metric
    outputfilename=sprintf('avgMod_AD_%s_%s',metrics,freq{f});
    range = [];
    range.minV = rangeVal.(freq{f}).AD(1); range.maxV = rangeVal.(freq{f}).AD(2);
    range.eps = 0; range.type = 2;
    [ctk]=[(othercolor('Mdarkrainbow',1000))]; CT = ctk;
    BNVplot(avg_AD_data,range,plttype,CT,odir,atlasname,outputfilename,adjustCM) 

    % plotting AD vs. Control
    %saving_only_contrast
    maxVal=max(max(survived_tval));
    minVal=min(min(survived_tval));
    if (maxVal==0 && minVal==0)
        disp('Not Significant')
    else
        outputfilename=sprintf('avgMod_AD_CONT_%s_%s_%s',metrics,freq{f},num2str(['GLMfit_FDR_p' afterPoint]));
        range = [];
        range.minV = rangeVal.(freq{f}).ADCONT(1); range.maxV = rangeVal.(freq{f}).ADCONT(2);
        range.eps = 0; range.type = 2;
        [ctk]=[flipud(othercolor('RdBu10',100))];CT = ctk;
        BNVplot(survived_tval,range,plttype,CT,odir,atlasname,outputfilename,adjustCM) 
    end  
end
outname=sprintf('Sig_val_modular_%s_%s',metrics,num2str(['GLMfit_FDR_p' afterPoint]));
filename_ = fullfile(outpath, outname);
save(filename_, '-struct','sig_vals'); 





