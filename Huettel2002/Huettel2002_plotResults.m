% Script to plot the results of simulation of Huettel et al (2002) data.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% Clear the place
clear; close('all');

% Load data (select data to plot)
% load Huettel2002_Simulation_Results
load Huettel2002_SimulationHMM_Results

figure(1); clf;
set(gcf, 'Color', ones(1,3));

for iObs = 1:3

    % PLOT RESULTS FOR REPETITIONS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subplot(4,2,3+(iObs-1)*2);
    errorbar((1:maxrep) - 0.1, mrep_viol(iObs,:), srep_viol(iObs,:), ...
        '.-', 'Color', [0 0 0],     'LineWidth', 2, 'MarkerSize', 25); 
    hold('on');
    errorbar((1:maxrep) + 0.1, mrep_cont(iObs,:), srep_cont(iObs,:), ...
        '.-', 'Color', 0.5*[1 1 1], 'LineWidth', 2, 'MarkerSize', 25);
    hold on
    errorbar((1:maxrep) - 0.1, mrep_viol(iObs,:), srep_viol(iObs,:), ...
        'Color', zeros(1,3), 'LineWidth', 2);
    xlabel('streak length'); ylabel({ObsName{iObs}; 'surprise'})
    ylim([0, 2.5]);
    set(gca, 'XTick', 1:8); box off
    xlim([0 maxrep+1])
    
    % PLOT RESULTS FOR ALTERNATIONS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subplot(4,2,4+(iObs-1)*2);
    errorbar((2:maxrep) - 0.1, malt_viol(iObs, 2:end), salt_viol(iObs, 2:end), ...
        '.-', 'Color', [0 0 0],     'LineWidth', 2, 'MarkerSize', 25); 
    hold('on');
    errorbar((2:maxrep) + 0.1, malt_cont(iObs, 2:end), salt_cont(iObs, 2:end), ...
        '.-', 'Color', 0.5*[1 1 1], 'LineWidth', 2, 'MarkerSize', 25);
    xlabel('streak length'); ylabel('surprise');
    ylim([0, 2.5]);
    set(gca, 'XTick', 1:8); box off
    xlim([0 maxrep+1])
end

% Plot original data by Huettel et al 2002.
load Huettel2002data.mat

subplot(4,2,1);
errorbar((1:8) - 0.1, rep_seq_viol_mean, rep_seq_viol_sem, '.-', 'Color', ...
    [0 0 0],     'LineWidth', 2, 'MarkerSize', 25); 
hold('on');
errorbar((1:8) + 0.1, rep_seq_cont_mean, rep_seq_cont_sem, '.-', 'Color', ...
    0.5*[1 1 1], 'LineWidth', 2, 'MarkerSize', 25);
legend({'violation', 'continuation'}, 'Location', 'NorthWest');
xlabel('streak length'); ylabel('RT (ms)'); title('REPETITIONS');
ylim([350 550]); 
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
ylim([350 550]); 
set(gca, 'XTick', 1:8); box off
xlim([0 8+1])
