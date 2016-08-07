% Script to plot the results of the simulations of Cho et al. (2002) data
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

clear all; close all;

% load data by Cho et al
original_data = load('Cho_et_al_data.mat');

% load simulation results (select which data to plot)
load Simulation_Cho_Results
% load Simulation_HMM_Cho_Results

% plot
figure(1); set(gcf, 'Color', [1 1 1]); set(gcf, 'Name', 'Model for Cho et al')

for iObs = 1:3
    
    % get data for this observer
    if iObs == 1; IOprediction_Cho_order = Prediction_TransitionProbs; end
    if iObs == 2; IOprediction_Cho_order = Prediction_StimProb; end
    if iObs == 3; IOprediction_Cho_order = Prediction_AltProb; end
    
    % plot the result
    subplot(2,2,1+iObs)
    plot(1:16, IOprediction_Cho_order, 'o-')
    xlim([0 17]); ylim([0.5 1.8]); grid on; box off       % for fixed belief model with leaky integration
    % xlim([0 17]); ylim([0.5 1.6]); grid on; box off           % for dynamic belief model
    ylabel({'surprise';'last observation'})
    if iObs == 1; title({'Inference based on';'transition probabilities'}); end
    if iObs == 2; title({'Inference based on';'stimulus frequency'}); end
    if iObs == 3; title({'Inference based on';'alternation frequency'}); end
    set(gca, 'XTick', 1:16, 'XTickLabel', patterns_Cho)
    set(gca, 'XTickLabelRotation', 90)
end

% plot data by Cho et al
subplot(2,2,1)
plot(1:16, original_data.RT, 'o-')
xlim([0 17]); ylim([290 430]); grid on; box off
ylabel('RT (ms)')
title('data by Cho et al 2002')
set(gca, 'XTick', 1:16, 'XTickLabel', patterns_Cho)
set(gca, 'XTickLabelRotation', 90)
