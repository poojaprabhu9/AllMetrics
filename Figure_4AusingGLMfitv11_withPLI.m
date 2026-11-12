% GLM fit for getting common regions across relative
% power+AECpair+Imcoh+PLI ==> peak_freq as covariate + age + gender +
% freq_labels * group_labels, a interaction term
%peak frequency used here is computed from 2 to 12 Hz
% plotting this for each metrics for mulpli correction p value=0.05/40*3*4
% [i.e. n_reg *n_freqh * n_metrics] and the only which survives after multiple comparison is used for common
% regions identification
clear all
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
figpath='C:/Users/pprabhu/OneDrive - UCSF/tresting/datafiles_PSD_dPTE/outputs_dell/Region_level_stats/Figures/commonReg_peakFreqLF2to12/';
datapath='C:/Users/pprabhu/OneDrive - UCSF/tresting/datafiles_PSD_dPTE/';
%%

metrics={'relative_pow','imcoh', 'aec_pair','pli'};
freq={'deltatheta','alpha','beta'};
p_multcmp=0.05/(n_reg*length(freq)*length(metrics));
% p_multcmp=0.05/(n_reg*length(freq));
p_thresh = 0.01;
%% Plotting parameters
plttype = 'full';
adjustCM = 1;
atlasname='BNmod';
range.deltatheta=[-50 50];

%% picking datafiles
AD=load(fullfile(datapath,'all_metrics_AD'));
CONT=load(fullfile(datapath,'all_metrics_CONT'));

pli_AD=load(fullfile(datapath,'pli_metrics_AD'));
pli_CONT=load(fullfile(datapath,'pli_metrics_CONT'));
AD.pli=pli_AD.pli;
CONT.pli=pli_CONT.pli;
peak_freq_AD_all=load(fullfile(datapath,'peak_freq_metrics_AD'));
peak_freq_CONT_all=load(fullfile(datapath,'peak_freq_metrics_CONT'));
peak_freq_LF1_AD= peak_freq_AD_all.modular.LF2to12;
% peak_freq_LF2_AD= peak_freq_AD_all.modular.LF7to12;
peak_freq_LF1_CONT= peak_freq_CONT_all.modular.LF2to12;
% peak_freq_LF2_CONT= peak_freq_CONT_all.modular.LF7to12;

peak_freq_AD= peak_freq_LF1_AD;
peak_freq_CONT= peak_freq_LF1_CONT;

%% dataSheets for demographic details
ADsheet = readtable(fullfile(datapath,"RadID_77AD_asPer_processed_data.xlsx"));
CONTsheet = readtable(fullfile(datapath,"RadID_90CONT_asPer_processed_data.xlsx"));

AD_age = ADsheet.Age_at_MEG;    
numericGender = cellfun(@(x) strcmp(x, 'Female'), ADsheet.Gender);  %female-1, male 0
AD_gender = double(numericGender);

CONT_age = CONTsheet.Age_at_MEG;    
numericGender = cellfun(@(x) strcmp(x, 'Female'), CONTsheet.Gender); 
CONT_gender = double(numericGender);

% group labels (4 for CONT and 5 for AD)
group_labelsTMP = [4 * ones(length(CONT_age), 1); 5 * ones(length(AD_age), 1)];
group_labels = [group_labelsTMP; group_labelsTMP; group_labelsTMP];

freq_labels = [7 * ones(length(group_labelsTMP), 1); 8 * ones(length(group_labelsTMP), 1); 9 * ones(length(group_labelsTMP), 1)]; % deltatheta, alpha, beta


age = [cat(1,CONT_age, AD_age); cat(1,CONT_age, AD_age); cat(1,CONT_age, AD_age)]; % three freq band
gender = [cat(1,CONT_gender, AD_gender); cat(1,CONT_gender, AD_gender); cat(1,CONT_gender, AD_gender)];

peak_freq_all = [cat(1,peak_freq_CONT, peak_freq_AD); cat(1,peak_freq_CONT, peak_freq_AD); cat(1,peak_freq_CONT, peak_freq_AD)];

%%
% results to store fitglm output
results = struct();
fc_data = [];
for f=1:length(freq)
    %%%%%%%%%%% AD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AD_rel_pow = AD.relative_pow.modular.(freq{f});
    AD_imcoh = squeeze(mean(AD.imcoh.modular.(freq{f}),2));
    AD_aec = squeeze(mean(AD.aec_pair.modular.(freq{f}),2));
    AD_pli = squeeze(mean(AD.pli.modular.(freq{f}),2));
    fc_data_AD = zeros(length(AD_rel_pow),n_reg,length(metrics) ); 
    for s=1:length(AD_rel_pow)
        for r=1:n_reg
            fc_data_AD(s,r,1)=AD_rel_pow(s,r);
            fc_data_AD(s,r,2)=AD_imcoh(s,r);
            fc_data_AD(s,r,3)=AD_aec(s,r); 
            fc_data_AD(s,r,4)=AD_pli(s,r);
        end
    end

    %%%%%%%%%%% Controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CONT_rel_pow = CONT.relative_pow.modular.(freq{f});
    CONT_imcoh = squeeze(mean(CONT.imcoh.modular.(freq{f}),2));
    CONT_aec = squeeze(mean(CONT.aec_pair.modular.(freq{f}),2));
    CONT_pli = squeeze(mean(CONT.pli.modular.(freq{f}),2));
    fc_data_CONT = zeros(length(CONT_rel_pow),n_reg,length(metrics)); 
    for s=1:length(CONT_rel_pow)
        for r=1:n_reg
            fc_data_CONT(s,r,1)=CONT_rel_pow(s,r);
            fc_data_CONT(s,r,2)=CONT_imcoh(s,r);
            fc_data_CONT(s,r,3)=CONT_aec(s,r); 
            fc_data_CONT(s,r,4)=CONT_pli(s,r);
        end
    end

     %%%%%%%%%%%% final functional connectivity data %%%%%%%%%%%%%%%
    fc_data_per_band = cat(1, fc_data_CONT, fc_data_AD);
    
    fc_data = cat(1,fc_data, fc_data_per_band);
    
end

for reg = 1:n_reg
        for m = 1:length(metrics)

           % dependent (or response) variable
            y=fc_data(:,reg,m);

            % peak_freq covariate
            peak_freq = peak_freq_all(:,reg);

            condition_labels = (m+1) * ones(length(fc_data),1);
    
            % covariates + dependent (or response) variable in table
            T = table(age, gender, freq_labels, group_labels, condition_labels, peak_freq, y);

            T.group_labels     = categorical(T.group_labels);
            T.freq_labels      = categorical(T.freq_labels); % enforce order
            T.condition_labels = categorical(T.condition_labels);
    
            % GLM analysis (assuming normal distribution)
            model_fit = fitglm(T, 'y ~ age + gender + peak_freq +  group_labels + freq_labels + group_labels * freq_labels +  condition_labels'); % Fit linear model with all covariates
    
            % store p-values and t-values
            results(reg,m).Coefficients = model_fit.Coefficients; 
        end
    end

%% Visualisation %%%%%%

clear significant_fitglmStats

tagname = {'group_labels_5','freq_labels_8:group_labels_5','freq_labels_9:group_labels_5'};
for f = 1:length(freq)
    significant_regions_perGroup= zeros(1,n_reg);
    significant_tvals_perGroup = zeros(1,n_reg);
    pval_perGroup_perMetric = ones(length(metrics),n_reg); % dont use zeros because 0 is less than 0.01 and even if its not significant it will show significant
    tval_perGroup_perMetric = zeros(length(metrics),n_reg);
    significant_pvals = ones(length(metrics),n_reg);
    for m =1:length(metrics)
        for reg = 1 : n_reg  
            pval_perGroup_perMetric(m,reg)=results(reg,m).Coefficients.pValue(tagname{f}); 
            tval_perGroup_perMetric(m,reg)=results(reg,m).Coefficients.tStat(tagname{f});          

       end
            significant_tvals = tval_perGroup_perMetric(m,:) .* (pval_perGroup_perMetric(m,:)<p_multcmp);
            significant_pvals(m,:) = pval_perGroup_perMetric(m,:) .* (pval_perGroup_perMetric(m,:)<p_multcmp);
            significant_pvals(m,(significant_pvals(m,:)==0)) = 1;
%             [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pval_perGroup_perMetric(m,:),p_thresh,'pdep','yes');
%             corrected_pval_perGroup_perMetric=adj_p.*h;
%             survived_tval_perGroup_perMetric=tval_perGroup_perMetric(m,:).*h;
%             significant_tvals = survived_tval_perGroup_perMetric;
%             clear sig_reg_tval
            sig_reg_tval = [significant_tvals(1,1:20) 0 0 0 0 significant_tvals(1,21:40) 0 0 0 0];
                  
            odir=figpath;
        
            range = [];
            rangeVal= max(abs(min(sig_reg_tval)),abs(max(sig_reg_tval)));
            range.minV = -rangeVal; range.maxV = rangeVal;
            range.eps = 0; range.type = 2;
            outputfilename=sprintf('common_BNmod_across_%s_AD_CONT_range%s_combinedFreq_%s_GLMfitAgeGenderPeakFreq2to12v11',metrics{m},num2str(range.maxV),freq{f});

            
            [ctk]=[flipud(othercolor('RdBu10',100))];CT = ctk;
            BNVplot(sig_reg_tval,range,plttype,CT,odir,atlasname,outputfilename,adjustCM)
    end
    for reg = 1:n_reg
%         pval_perGroup = pval_perGroup_perMetric(:,reg);
        pval_perGroup = significant_pvals(:,reg); % pvalues that survived after p_multcmp
        if all(pval_perGroup<p_thresh)
            significant_regions_perGroup(1,reg)=1;
            significant_tvals_perGroup(1,reg)= sum(tval_perGroup_perMetric(:,reg),'all');
        end
    end
    significant_fitglmStats.(freq{f}).CommonRegions=significant_regions_perGroup;
    significant_fitglmStats.(freq{f}).tvals=significant_tvals_perGroup;

    clear sig_reg_tval
    sig_reg_tval = [significant_tvals_perGroup(1,1:20) 0 0 0 0 significant_tvals_perGroup(1,21:40) 0 0 0 0];
    significant_fitglmStats.(freq{f}).tvals48Reg = sig_reg_tval;

    odir=figpath;
    outputfilename=sprintf('common_BNmod_across_4Metrics_AD_CONT_range%s_combinedFreq_%s_GLMfitAgeGenderPeakFreq2to12v11',num2str(range.maxV),freq{f});

    range = [];
    rangeVal= max(abs(min(sig_reg_tval)),abs(max(sig_reg_tval)));
    range.minV = -rangeVal; range.maxV = rangeVal;
    range.eps = 0; range.type = 2;
    
    [ctk]=[flipud(othercolor('RdBu10',100))];CT = ctk;
    BNVplot(sig_reg_tval,range,plttype,CT,odir,atlasname,outputfilename,adjustCM) 
end
% save
foutname=sprintf('T_val_%s_AD_modular_4Metrics_combinedFreq_GLMfitAgeGenderPeakFreq2to12v11',num2str(p_thresh));
filename_ = fullfile(outpath, foutname);
save(filename_,'significant_fitglmStats'); 