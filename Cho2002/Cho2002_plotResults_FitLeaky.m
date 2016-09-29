% Script to plot the results of the simulations of Cho et al. (2002) data
% The plot shows the results for each observer, and for the best-fitting
% leak value.
%
% Copyright 2016 Florent Meyniel & Maxime Maheu 

clear all; close all;

% load data by Cho et al
original_data = load('Cho2002data.mat');

% load simulation results
load Cho2002_Simulation_Results_FitLeaky

% prepare plot
figure(1); set(gcf, 'Color', [1 1 1]); set(gcf, 'Name', 'Model for Cho et al')

% plot data by Cho et al
figure(1)
subplot(2,2,1)
plot(1:16, original_data.RT, 'o-')
RT_lim = [290 430];
xlim([0 17]); ylim(RT_lim); grid on; box off
ylabel('RT (ms)')
title('data by Cho et al 2002')
set(gca, 'XTick', 1:16, 'XTickLabel', patterns_Cho)
set(gca, 'XTickLabelRotation', 90)

% plot ideal observer prediction
allR2 = nan(3, numel(grid_leak));
allBeta = nan(3, numel(grid_leak), 2);
for iObs = 1:3
    
    % get data for this observer
    if iObs == 1; LearnedStat = 'transition probabilities'; IOprediction_Cho_order = Prediction_TransitionProbs; end
    if iObs == 2; LearnedStat = 'stimulus frequency' ; IOprediction_Cho_order = Prediction_StimProb; end
    if iObs == 3; LearnedStat = 'alternation frequency' ; IOprediction_Cho_order = Prediction_AltProb; end
    
    % find the best fitting parameter value
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for k = 1:numel(grid_leak)
        out = regstats(original_data.RT, IOprediction_Cho_order(k,:)', 'linear', {'beta','rsquare'});
        allR2(iObs, k) = out.rsquare;
        allBeta(iObs, k, :) = out.beta;
    end
    [~, best_k] = max(allR2(iObs,:));
    fprintf('\n %s best-fitting k=%d, R2: %3.2f',LearnedStat, grid_leak(best_k), allR2(iObs, best_k))
    
    % scale axis to match RTs
    Surprise_lim = (RT_lim - allBeta(iObs, best_k,1)) / allBeta(iObs, best_k,2);
    
    % plot the result
    subplot(2,2,1+iObs)
    plot(1:16, IOprediction_Cho_order(best_k,:), 'o-')
    xlim([0 17]); ylim(Surprise_lim); grid on; box off       % for fixed belief model with leaky integration
    ylabel({'surprise';'last observation'})
    title({'Inference based on';LearnedStat})
    set(gca, 'XTick', 1:16, 'XTickLabel', patterns_Cho)
    set(gca, 'XTickLabelRotation', 90)
end

