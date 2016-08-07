% Script to plot the results of the simulations of Falk (1975).
% The simulation corresponds to "dynamic belief models".
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

clear

% LOAD SIMULATION RESULTS
load Comparison_NonStationarity_data

% Get the prior probability on changes
pJump_str = cell(1, numel(pJump));
for k = 1:numel(pJump)
    pJump_str{k} = ['1/', num2str(1/pJump(k))];
end

% PLOT RESULTS
figure(2); set(gcf, 'Color', [1 1 1]); clf
set(gcf, 'Name', 'Ideal Observers with jump, different prior believes on volatility')
col = winter(numel(pJump));

subplot(2,2,1)
Falkdata = load('Falk1975data.mat');
plot(Falkdata.prob_of_alternation, Falkdata.perceived_randomness, 'ko-');
grid on; box off
xlabel('Prob. of Alternation')
ylabel('perceived randomness')
title('data by Falk 1975')

subplot(2,2,2)
hold on
for iP = 1:numel(pJump)
    plot(pAlt, m_avg_entropy_trans(:,iP), 'color', col(iP,:), 'LineWidth', 2)
end
legend(pJump_str, 'location', 'South')
ylim([0.5 1])
grid on
xlabel('Prob. of Alternation')
ylabel('average posterior entropy')
title('When learning transition probs.')

subplot(2,2,3)
hold on
for iP = 1:numel(pJump)
    plot(pAlt, m_avg_entropy_pAlt(:,iP), 'color', col(iP,:), 'LineWidth', 2)
end
legend(pJump_str, 'location', 'South')
ylim([0.5 1])
grid on
xlabel('Prob. of Alternation')
ylabel('average posterior entropy')
title('When learning p(Alt)')

subplot(2,2,4)
hold on

for iP = 1:numel(pJump)
    plot(pAlt, m_avg_entropy_freq(:,iP), 'color', col(iP,:), 'LineWidth', 2)
end
legend(pJump_str, 'location', 'South')
ylim([0.5 1])
grid on
xlabel('Prob. of Alternation')
ylabel('average posterior entropy')
title('When learning p(stim)')


