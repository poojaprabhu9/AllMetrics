%correlation between different metrics and MMSE
% used on common regions tval for each metrics
% excluded dPTE
clear all
close all
outpath='C:/Users/pprabhu/OneDrive - UCSF/tresting/datafiles_PSD_dPTE/outputs_dell/Region_level_stats/sig_values/';
% data=readtable(fullfile(outpath,'datasheet_after_reordering_AD.xlsx')); % close the excel file in the background
data=readtable(fullfile(outpath,'datasheet_after_reordering_AD_common_regions_glmfitAgeGender.xlsx')); % close the excel file in the background
stats_criteria= 'FDR001'; 
% tval_common =load(fullfile(outpath,['T_val_Common_regions_modular_3metrics_' stats_criteria]));

%%
metrics={'relative_pow','aec_pair','imcoh'}; % 'relative_pow','imcoh', 'pte','aec_pair','peak_freq'
freq={'deltatheta','alpha','beta'};
measure = {'MMSE','CDR_BOX'};
colorVal = [195/255 64/255 127/255; 0/255 79/255 81/255]; % for each measure
colorVal = [228/255 87/255 46/255; 0/255 79/255 81/255]; % for each measure


bonferr_p_thresh = 0.05/3; %across 3 metrics 

metricsPlot1 = getPlotHandles(1,2,[0.05 0.6 0.28 0.3],0.015,0.06);
ytickVal1 = [10 50; 0 0.15; 0.02 0.05];

metricsPlot2 = getPlotHandles(1,2,[0.38 0.6 0.28 0.3],0.015,0.06);
ytickVal2 = [10 35; 0 0.07; 0.03 0.06];

metricsPlot3 = getPlotHandles(1,2,[0.71 0.6 0.28 0.3],0.015,0.06);
ytickVal3 = [10 50; 0 0.06; 0.01 0.03];

AD_age = data.Age_at_MEG;    
numericGender = cellfun(@(x) strcmp(x, 'Female'), data.Gender);  %female-1, male 0
AD_gender = double(numericGender);

age=[AD_age; AD_age; AD_age];
gender=[AD_gender; AD_gender; AD_gender];
gender = categorical(gender);

metrics_label=[ones(length(AD_age),1); 2*ones(length(AD_age),1); 3*ones(length(AD_age),1)];
metrics_label = categorical(metrics_label);

% SubjectID
subject = [1:length(AD_age) 1:length(AD_age) 1:length(AD_age)]'; % copied three times for three metrics

for f=1:length(freq)
    if f==1
        metricsPlot = metricsPlot1;
        ytickVal = ytickVal1;

    elseif f==2
        metricsPlot = metricsPlot2;
        ytickVal = ytickVal2;

    elseif f==3
        metricsPlot = metricsPlot3;
        ytickVal = ytickVal3;

    end
%     colorVal = [0.5 0.5 0.5];
%     colorVal = [0 0 0];
%     for m=1:length(metrics)
        fieldname1=sprintf('%s_common_reg_%s',metrics{1},freq{f});
        fieldname2=sprintf('%s_common_reg_%s',metrics{2},freq{f});
        fieldname3=sprintf('%s_common_reg_%s',metrics{3},freq{f});
        metricsVals=[data.(fieldname1); data.(fieldname2); data.(fieldname3)];
        
        for ms =1:length(measure)
            y=[data.(measure{ms}); data.(measure{ms}); data.(measure{ms})];

            T = table(y, metricsVals, metrics_label, age, gender, subject);

%             % GLM analysis (assuming normal distribution)
%             model_fit = fitglm(T, 'y ~ metricsVals + metrics_label + age + gender '); % Fit linear model with all covariates
            
          
            % Linear mixed model formula
            % MMSE ~ metric Values + categorical(metrics label)  + age + categorical(gender) + (1|subject)
            formula = 'y ~ metricsVals + metrics_label + gender + age + (1|subject)';
            
            % Fit the linear mixed model
            model_fit = fitlme(T, formula);

            % predicted Measure (MMSE or CDR_BOX) values 
            predicted_y = predict(model_fit);

            % add back to table
            T.predicted_measure = predicted_y;

            % Scatter plot of actual vs predicted MMSE
            scatter(metricsPlot(1,ms),T.y, T.predicted_measure, 10, 'filled', 'MarkerEdgeColor', [colorVal(ms,1) colorVal(ms,2) colorVal(ms,3)], 'MarkerFaceColor', [colorVal(ms,1) colorVal(ms,2) colorVal(ms,3)],'Marker','o');
            xlabel(metricsPlot(1,ms), ['Actual ' measure{ms}]);
            ylabel(metricsPlot(1,ms), ['Predicted ' measure{ms}]);
            hold(metricsPlot(1,ms),"on")
            
            % Linear regression model to the data for the fitted line
            model_lineFit = fitlm(T.y, T.predicted_measure);
            
            % Get the predicted values and confidence intervals
            [ypred, delta] = predict( model_lineFit, T.y);
            
            % Plot the fitted line
            plot(metricsPlot(1,ms), T.y, ypred, 'k-', 'LineWidth', 2); 
            
            % Plot the confidence interval
            plot(T.y, ypred + delta, 'k--', 'LineWidth', 1);  % Upper bound of CI
            plot(T.y, ypred - delta, 'k--', 'LineWidth', 1);  % Lower bound of CI
            
            set(metricsPlot(1,ms),'XGrid','off','YGrid','on')
%             legend(metricsPlot(1,ms),'Actual vs Predicted', 'Fitted Line', '95% Confidence Interval', 'Location', 'best');

        end
end
