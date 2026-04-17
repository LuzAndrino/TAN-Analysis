
%Creating the vector - pulled from excel sheet for now 
% - data from Mar 2,2026
% young_synapses = [28.70, 23.62, 22.98, 25.92, 22.78, 24.87, 23.51, 21.33, 27.73, 22.29, 19.28];
% young_efr      = [1.618, 1.456, 1.755, 3.504, 1.814, 0.8171, 1.442, 1.08, 2.045, 4.619, 8.845];
% 
% old_synapses = [22.08, 18.40, 21.67, 22.78, 22.64, 18.64, 23.27, 16.29, 20.15];
% old_efr      = [0.8281, 0.6846, 0.06367, 1.047, 0.7874, 0.9542, 0.7863, 0.1312, 0.4364];

%%
%Young
young_ID = 
young_synapses = 
young_efr =
%% make into useable tables
 young_table = table(young_synapses', young_efr', ...
    'VariableNames', {'Synapses', 'EFR_SL'});

  old_table = table(old_synapses', old_efr', ...
    'VariableNames', {'Synapses', 'EFR_SL'});

young_table(11, :) = []; %outlier detected at EFR 8.845 so I drop full row 

%% correlation stats 
[r_young, p_young] = corr(young_table.Synapses, young_table.EFR_SL, 'Type', 'Spearman');
[r_old,   p_old]   = corr(old_table.Synapses,   old_table.EFR_SL,   'Type', 'Spearman');

all_synapses = [young_table.Synapses; old_table.Synapses]; % combine data for input matrix for polypredci function
all_efr      = [young_table.EFR_SL;   old_table.EFR_SL];
[r_all, p_all] = corr(all_synapses, all_efr, 'Type', 'Spearman');

%% r an d p values reported 
fprintf('Young:    r = %.3f, p = %.4f\n', r_young, p_young)
fprintf('Old:      r = %.3f, p = %.4f\n', r_old,   p_old)
fprintf('Combined: r = %.3f, p = %.4f\n', r_all,   p_all)

%% Plotting- Correlation between synapses and EFR SL
% Plot dots by group
figure; hold on;

inputmatrix = [all_synapses'; all_efr'];

Scatterplot_polypredci(inputmatrix);
delete(findobj(gca, 'Marker', 'o'))

scatter(young_table.Synapses, young_table.EFR_SL, 60, [0.99 0.55 0.35], 'filled')  %orange
scatter(old_table.Synapses,   old_table.EFR_SL,   60, [0.42 0.68 0.84], 'filled')  %  blue

% One trendline through ALL data
% all_synapses = [young_synapses, old_synapses];
% all_efr      = [young_efr,      old_efr];
% 
% p = polyfit(all_synapses, all_efr, 1);  % linear fit

% x_line = linspace(min(all_synapses)-1, max(all_synapses)+1, 100);
% y_line = polyval(p, x_line);
% plot(x_line, y_line, 'k-', 'LineWidth', 1.5)

% Labels and formatting
xlabel('Synapse Count@3kHz', 'FontSize', 13)
ylabel('30dB EFR SL Amplitude (\muV)', 'FontSize', 13)
legend({'Trendline', 'Young', 'Old'}, 'Location', 'northwest')
set(gca, 'box', 'off', 'FontSize', 12)  % cleaner look

% function outliers = tukey_outliers(data, label)
%     Q1  = quantile(data, 0.25);
%     Q3  = quantile(data, 0.75);
%     IQR = Q3 - Q1;
%     lower = Q1 - 1.5 * IQR;
%     upper = Q3 + 1.5 * IQR;
% 
%     outliers = data(data < lower | data > upper);
% 
%     fprintf('%s -- Lower fence: %.3f | Upper fence: %.3f\n', label, lower, upper)
%     if isempty(outliers)
%         fprintf('%s -- No outliers detected\n\n', label)
%     else
%         fprintf('%s -- Outliers: ', label)
%         disp(outliers)
%     end
% end
% 
% % Run on each group
% tukey_outliers(young_efr, 'Young EFR') %Young EFR -- Outliers:     4.6190 (1,12)   8.8450 (1,13)
% tukey_outliers(old_efr,   'Old EFR')
% tukey_outliers(young_synapses, 'Young Synapses')
% tukey_outliers(old_synapses,   'Old Synapses')