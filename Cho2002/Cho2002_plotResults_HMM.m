% Script to plot the results of the simulations of Cho et al. (2002) data
% This code is for the "dynamic belief models" (using HMM). The prior
% volatility is not fitted, instead, we took the best-fitting value 
% obtained with Squires et al data.
%
% Copyright 2016 Florent Meyniel & Maxime Maheu 

clear all; close all;

% load data by Cho et al
original_data = load('Cho2002data.mat');

% load simulation results
load Cho2002_Simulation_Results_HMM

% plot
figure(1); set(gcf, 'Color', [1 1 1]); set(gcf, 'Name', 'Model for Cho et al')

for iObs = 1:3
    
    % get data for this observer
    if iObs == 1; LearnedStat = 'transition probabilities'; IOprediction_Cho_order = Prediction_TransitionProbs; end
    if iObs == 2; LearnedStat = 'stimulus frequency' ; IOprediction_Cho_order = Prediction_StimProb; end
    if iObs == 3; LearnedStat = 'alternation frequency' ; IOprediction_Cho_order = Prediction_AltProb; end
    
    % plot the result
    subplot(2,2,1+iObs)
    plot(1:16, IOprediction_Cho_order, 'o-')
    xlim([0 17]); ylim([0.5 1.6]); grid on; box off           % for dynamic belief model
    ylabel({'surprise';'last observation'})
    title({'Inference based on';LearnedStat})
    set(gca, 'XTick', 1:16, 'XTickLabel', patterns_Cho)
    set(gca, 'XTickLabelRotation', 90)
    
    % report goodness of the fit
    out = regstats(original_data.RT, IOprediction_Cho_order, 'linear', {'rsquare', 'yhat'});
    fprintf('\n %s R2: %3.2f',LearnedStat, out.rsquare)
end

% plot data by Cho et al
subplot(2,2,1)
plot(1:16, original_data.RT, 'o-')
xlim([0 17]); ylim([290 430]); grid on; box off
ylabel('RT (ms)')
title('data by Cho et al 2002')
set(gca, 'XTick', 1:16, 'XTickLabel', patterns_Cho)
set(gca, 'XTickLabelRotation', 90)
