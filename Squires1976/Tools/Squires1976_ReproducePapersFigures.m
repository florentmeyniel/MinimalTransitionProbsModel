% This script reproduce some of the figures in the Squires et al. (1976)
% and Kolossa et al. (2003)
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% Clear the place
close('all'); clear;

%% REPRODUCE THE FIGURE 8 OF KOLOSSA ET AL. (2003)
%  ===============================================

% Load data
[KolossaData, nPat, pAlist] = ...
    Squires1976_PrepareData('DataKolossaFIHN2013.mat');

% Reorder the blocks
[~,I] = ismember([0.7, 0.5, 0.3], pAlist);
KolossaData = KolossaData(:,:,I);
pAlist = pAlist(I);

% Figure options
figure('Name', 'Kolossa et al. (2003) - Figure 8', ...
    'Units', 'Normalized', 'Position', [.7 .25 .25 .5]);
col    = 'g';
stim   = [2, 1, 1];
letter = {'A', 'B'};

% For each block
for p = 1:numel(pAlist)
    
    % Plot the tree
    subplot(1,3,p);
    PlotTree(KolossaData(:,:,p), nPat, col, [], [], stim(p));
    
    % Plot options
    set(gca, 'Box', 'Off', 'XTickLabel', 3:-1:0);
    axis([0, 3, 1.25, 6]);
    xlabel('order', 'Interpreter', 'LaTeX', 'FontSize', 15);
    if p == 1, ylabel('P300 amplitude [$\mu$V]', 'Interpreter', 'LaTeX', 'FontSize', 15); end
    title(sprintf('$p$(%s) = %1.2f', letter{stim(p)}, pAlist(p)), ...
        'Interpreter', 'LaTeX', 'FontSize', 15);
    
    % For each element on the plot
    c = get(gca, 'Children');
    for i = 1:numel(c)
        
        % Change the color
        set(c(i), 'Color', col, 'LineWidth', 1);
        
        % Customize the labels
        if ~isempty(strfind(class(c(i)), 'Text'))
            set(c(i), 'String', ['\it{',lower(get(c(i), 'String')),'}'], ...
                'FontSize', 15, 'Interpreter', 'LaTeX', 'Color', 'k');
        
        % Customize the edges and the dots
        else
            set(c(i), 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'none');
        end
    end
end

%% REPRODUCE THE FIGURE 1 OF SQUIRES ET AL. (1976)
%  ===============================================

% Load data
[SquiresData, nPat, pAlist, PatLabels] = ...
    Squires1976_PrepareData('DataSquiresScience1976.mat');

% Figure options
figure('Name', 'Squires et al. (1976) - Figure 1', ...
    'Units', 'Normalized', 'Position', [.3 .25 .4 .5]);

% For each block
for p = 1:numel(pAlist)
    
    % Plot the tree
    subplot(1,numel(pAlist),p);
    PlotTree(SquiresData(:,:,p), nPat, 'k', [], [], 1, 6, 1);
    
    % Plot options
    set(gca, 'YTick', [-0.5, 0:5, 5.5], 'Box', 'Off', ...
        'XMinorTick', 'off', 'YMinorTick', 'on');
    set(gca, 'YTickLabel', cellfun(@(x) sprintf('%1.1f', x), ...
        num2cell(get(gca, 'YTick')), 'UniformOutput', false));
    ylim([-0.5, 5.5]);
    if p > 1, set(gca, 'YColor', 'w');
    else ylabel('Discriminant score'); end
    xlabel('Order', 'FontSize', 10);
    title({sprintf('p(A) = %1.1f', pAlist(p)),''}, 'FontSize', 10);
    
    % For each element on the plot
    c = get(gca, 'Children');
    for i = 1:numel(c)
        if ~isempty(strfind(class(c(i)), 'Line'))
            set(c(i), 'Marker', 'None');
        end
    end
end

%% REPRODUCE THE FIGURE 3 OF SQUIRES ET AL. (1976)
%  ===============================================

% Reorder p(A) so that they are in the same order as in the paper
[pAlist,I] = sort(pAlist);
SquiresData = SquiresData(:,:,I);

% Generate one sequence per block
nPat = size(SquiresData,1);
L = 200;
Seq = Squires1976_GenerateSequences(pAlist, L, 1, nPat);

% Prepare the figure
figure('Name', 'Squires et al. (1976) - Figure 3', ...
    'Units', 'Normalized', 'Position', [.05 .25 .25 .5]); hold('on');

% For each block (defined by a global probability)
for b = 1:numel(pAlist)
    
    % Get the P300 amplitude
    DiscriminantScore = SquiresData(end,:,b); % y-axis
    
    % Get the prediction of the Squires et al. (1976) model
    s = Seq(1,:,b);
    Expectancy = Squires1976_SquiresModel(s, pAlist(b), nPat);
    Expectancy = GetPatLL(nPat, s, Expectancy);
    Expectancy = Expectancy(end,:);
    
    % Plot the results
    if     b == 1
        plot(Expectancy, DiscriminantScore, '^', 'MarkerSize', 6, ...
            'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    elseif b == 2
        plot(Expectancy, DiscriminantScore, 'ko', 'MarkerSize', 8, 'LineWidth', 1);
        plot(Expectancy, DiscriminantScore, 'k.');
    elseif b == 3
        plot(Expectancy, DiscriminantScore, 's', 'MarkerSize', 7, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1);
        plot(Expectancy, DiscriminantScore, 's', 'MarkerSize', 12, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1);
        plot(Expectancy, DiscriminantScore, 'k.');
    end
end

% Plot options
xlabel('Expectancy');
ylabel('Discriminant score');
axis([0.1, 0.7, -1, 6]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
