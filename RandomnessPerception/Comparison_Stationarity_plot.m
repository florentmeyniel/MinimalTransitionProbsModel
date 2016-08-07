% Script to plot the results of the simulations of Falk (1975).
% The simulation corresponds to "fixed belief models" with leaky integration.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

clear

% LOAD SIMULATION RESULTS
load Comparison_Stationarity_data

% get length of sequences
memDecay_str = cell(1, numel(memDecay));
for k = 1:numel(memDecay)
    memDecay_str{k} = num2str(memDecay(k));
end

% PLOT RESULTS
figure(1); set(gcf, 'Color', [1 1 1]); clf
set(gcf, 'Name', 'Ideal Observer without jump, different sequence lengths')
col = winter(numel(memDecay));

subplot(2,2,1)
Falkdata = load('Falk1975data.mat');
plot(Falkdata.prob_of_alternation, Falkdata.perceived_randomness, 'ko-');
grid on; box off
xlabel('Prob. of Alternation')
ylabel('perceived randomness')
title('data by Falk 1975')

subplot(2,2,2)
hold on
for iD = 1:numel(memDecay)
    plot(pAlt, m_avg_entropy_trans(:,iD), 'color', col(iD,:), 'LineWidth', 2)
end
legend(memDecay_str, 'location', 'South')
ylim([0.5 1])
grid on
xlabel('Prob. of Alternation')
ylabel('average posterior entropy')
title('When learning transition probs.')

subplot(2,2,3)
hold on
for iD = 1:numel(memDecay)
    plot(pAlt, m_avg_entropy_pAlt(:,iD), 'color', col(iD,:), 'LineWidth', 2)
end
legend(memDecay_str, 'location', 'South')
ylim([0.5 1])
grid on
xlabel('Prob. of Alternation')
ylabel('average posterior entropy')
title('When learning p(Alt)')

subplot(2,2,4)
hold on
for iD = 1:numel(memDecay)
    plot(pAlt, m_avg_entropy_freq(:,iD), 'color', col(iD,:), 'LineWidth', 2)
end
legend(memDecay_str, 'location', 'South')
ylim([0.5 1])
grid on
xlabel('Prob. of Alternation')
ylabel('average posterior entropy')
title('When learning p(stim)')
