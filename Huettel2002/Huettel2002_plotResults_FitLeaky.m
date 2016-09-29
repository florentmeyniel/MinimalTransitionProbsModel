% Script to plot the results of simulation of Huettel et al (2002) data.
% The plot shows the results for each observer, and for the best-fitting
% leak value.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% Clear the place
clear; % close('all');

% Load data (select data to plot)
% load Huettel2002_Simulation_Results
% load Huettel2002_SimulationHMM_Results
load('Huettel2002_Simulation_Results_2016-9-23_17-36-59_29387.mat')


% set the stage
figure(1); clf;
set(gcf, 'Color', ones(1,3));

% PLOT ORIGINAL DATA BY HUETTEL ET AL 2002.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load Huettel2002data.mat
RT_lim = [350 550];

subplot(4,2,1);
errorbar((1:8) - 0.1, rep_seq_viol_mean, rep_seq_viol_sem, '.-', 'Color', ...
    [0 0 0],     'LineWidth', 2, 'MarkerSize', 25); 
hold('on');
errorbar((1:8) + 0.1, rep_seq_cont_mean, rep_seq_cont_sem, '.-', 'Color', ...
    0.5*[1 1 1], 'LineWidth', 2, 'MarkerSize', 25);
legend({'violation', 'continuation'}, 'Location', 'NorthWest');
xlabel('streak length'); ylabel('RT (ms)'); title('REPETITIONS');
ylim(RT_lim); 
set(gca, 'XTick', 1:8); box off
xlim([0 8+1])

subplot(4,2,2);
errorbar((1:8) - 0.1, alt_seq_viol_mean, alt_seq_viol_sem, '.-', 'Color', ...
    [0 0 0],     'LineWidth', 2, 'MarkerSize', 25); 
hold('on');
errorbar((1:8) + 0.1, alt_seq_cont_mean, alt_seq_cont_sem, '.-', 'Color', ...
    0.5*[1 1 1], 'LineWidth', 2, 'MarkerSize', 25);
legend({'violation', 'continuation'}, 'Location', 'NorthWest');
xlabel('streak length'); ylabel('RT (ms)'); title('ALTERNATIONS');
ylim(RT_lim); 
set(gca, 'XTick', 1:8); box off
xlim([0 8+1])

% PLOT DATA FOR THE IDEAL OBSERVER
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
allR2 = nan(3, numel(grid_leak));
allBeta = nan(3, numel(grid_leak), 2);
for iObs = 1:3
    
    % FIND THE BEST FITTING PARAMETER VALUE
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for k = 1:25
        % Report the goodness of fit
        data = [rep_seq_viol_mean; rep_seq_cont_mean; alt_seq_viol_mean(2:end); alt_seq_cont_mean(2:end)];
        pred = [squeeze(mrep_viol(iObs,:,k))'; squeeze(mrep_cont(iObs,:,k))'; ...
            squeeze(malt_viol(iObs,2:end,k))'; squeeze(malt_cont(iObs,2:end,k))'];
        
        out = regstats(data, pred, 'linear', {'beta','rsquare'});
        allR2(iObs, k) = out.rsquare;
        allBeta(iObs, k, :) = out.beta;
    end
    
    % find the index of the best-fitting parameter
    [~, best_k] = max(allR2(iObs,:));
    
    % Report best fitting value and goodness-of-fit
    fprintf('\n %s, best-fitting k=%d, R2: %3.2f',ObsName{iObs}, grid_leak(best_k), allR2(iObs, best_k))
    
    % scale axis to match RTs
    Surprise_lim = (RT_lim - allBeta(iObs, best_k,1)) / allBeta(iObs, best_k,2);
    
    % PLOT RESULTS FOR REPETITIONS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subplot(4,2,3+(iObs-1)*2);
    errorbar((1:maxrep) - 0.1, mrep_viol(iObs,:,best_k), srep_viol(iObs,:,best_k), ...
        '.-', 'Color', [0 0 0],     'LineWidth', 2, 'MarkerSize', 25); 
    hold('on');
    errorbar((1:maxrep) + 0.1, mrep_cont(iObs,:,best_k), srep_cont(iObs,:,best_k), ...
        '.-', 'Color', 0.5*[1 1 1], 'LineWidth', 2, 'MarkerSize', 25);
    hold on
    errorbar((1:maxrep) - 0.1, mrep_viol(iObs,:,best_k), srep_viol(iObs,:,best_k), ...
        'Color', zeros(1,3), 'LineWidth', 2);
    xlabel('streak length'); ylabel({ObsName{iObs}; 'surprise'})
    ylim(Surprise_lim);
    set(gca, 'XTick', 1:8); box off
    xlim([0 maxrep+1])
    
    % PLOT RESULTS FOR ALTERNATIONS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subplot(4,2,4+(iObs-1)*2);
    errorbar((2:maxrep) - 0.1, malt_viol(iObs, 2:end,best_k), salt_viol(iObs, 2:end,best_k), ...
        '.-', 'Color', [0 0 0],     'LineWidth', 2, 'MarkerSize', 25); 
    hold('on');
    errorbar((2:maxrep) + 0.1, malt_cont(iObs, 2:end,best_k), salt_cont(iObs, 2:end,best_k), ...
        '.-', 'Color', 0.5*[1 1 1], 'LineWidth', 2, 'MarkerSize', 25);
    xlabel('streak length'); ylabel('surprise');
    ylim(Surprise_lim);
    set(gca, 'XTick', 1:8); box off
    xlim([0 maxrep+1])
end


