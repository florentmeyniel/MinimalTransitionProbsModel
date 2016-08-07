% Plot the best fitting model.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 

%% INITIALIZATION
%  ==============

% Clear the place
clear; close('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE SOME OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% On which dataset
d = 'SquiresScience1976';
%d = 'KolossaFIHN2013';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the functions
addpath(genpath('Squires1976'));
try cd('Squires1976'); catch, end

%% LOAD DATA
%  =========

% Load results from the best models only
load(sprintf('Squires1976_BestParameters_%s.mat', d))
stat = Info.Stat;
nStat = numel(stat);
intg = Info.Intg;
nIntg = numel(intg);

% Load the data
[Data, nPat, pX, PatLabs] = Squires1976_PrepareData(sprintf('Data%s.mat', d));

% Reorder the blocks
[~,I] = sort(pX, 2, 'descend');
Data = Data(:,:,I);
pX = pX(I);
nP = numel(pX);

% Get patterns' expected frequency
PatProb = NaN(size(Data));
for s = 1:nP
    PatProb(:,:,s) = GetPatProb(nPat, pX(s));
end

% Font size and line width for the plots
fs = 12;
lw = 1.5;

% Colors to use
scol = [255 000 000;         % X
        000 000 255] ./ 255; % Y
ocol = [scol(1,:); mean(scol,1); scol(2,:)];
col = [(ocol(1,:) .* pX(1)); ...
        ocol(2,:); ...
       (ocol(end,:) .* pX(1))];
g   = repmat(128,1,3)./255;

% Define the name of the dataset
if     strcmp(d, 'SquiresScience1976'), dataset = {'Dataset from Squires et al. (1976)', ''};
elseif strcmp(d, 'KolossaFIHN2013'),    dataset = {'Dataset from Kolossa et al. (2013)', ''};
end

%% PLOT THE QUALITY OF FIT FOR THE BEST PARAMETERS
%  ===============================================

% Display the MSE
Var = 'BIC';
figure('Units', 'Normalized', 'Position', [0, 0.77, 0.25, 0.38]);
subplot(2,1,1);
h = bar([Best.MSE(:); Desc.MSE'], 'k', 'LineWidth', lw);

% Customize the axes
xlim([0,numel(Comp.(Var).All.weights)+1]);
set(gca, 'Box', 'Off', 'Layer', 'Bottom', 'LineWidth', lw, 'FontSize', fs);
ylabel('Mean squared error', 'FontSize', fs);

%
subplot(2,1,2);
bar(Comp.(Var).All.weights, 'k');
for i = 1:numel(Comp.(Var).All.weights)
    text(i, Comp.(Var).All.weights(i)+0.1, sprintf('%2.0f%%', Comp.(Var).All.weights(i)*100));
end
axis([0,numel(Comp.(Var).All.weights)+1,0,1]);
set(gca, 'Box', 'Off');
ylabel(sprintf('%s weights', Var));

%% DISPLAY THE RESULTS OF THE MODEL COMPARISON
%  ===========================================

%
figure('Position', [0.28 0.55 0.4 0.38]);
subplot(2,2,2);
barh(Comp.(Var).BayesStat.weights, 0.6, 'k', 'LineWidth', lw); hold('on');
plot(repmat(1/numel(Info.Stat),1,2), [0,numel(Info.Stat)+0.5], '-', 'Color', repmat(128/255,1,3), 'LineWidth', lw);
for i = 1:numel(Info.Stat)
    text(Comp.(Var).BayesStat.weights(i)+0.1, i, sprintf('%2.0f%%', Comp.(Var).BayesStat.weights(i)*100));
end
axis([0,1,0.5,numel(Info.Stat)+0.5]); axis('square');
set(gca, 'YTickLabel', {}, 'XAxisLocation', 'Top', 'Box', 'Off', 'Layer', 'Top', 'LineWidth', lw, 'Ydir', 'Reverse');
title({sprintf('%s family weights', Var),''});

%
subplot(2,2,4);
barh(Comp.(Var).BayesIntg.weights, 0.8, 'k', 'LineWidth', lw); hold('on');
plot(repmat(1/numel(Info.Intg),1,2), [0,numel(Info.Intg)+0.5], '-', 'Color', repmat(128/255,1,3), 'LineWidth', lw);
for i = 1:numel(Info.Intg)
    text(Comp.(Var).BayesIntg.weights(i)+0.1, i, sprintf('%2.0f%%', Comp.(Var).BayesIntg.weights(i)*100));
end
axis([0,1,0.5,numel(Info.Intg)+0.5]); axis('square');
set(gca, 'YTickLabel', {}, 'XAxisLocation', 'Top', 'Box', 'Off', 'Layer', 'Top', 'LineWidth', lw, 'Ydir', 'Reverse');

%
subplot(2,2,1);
imagesc(Comp.(Var).BayesStat.normevidence); hold('on'); colormap(flipud(autumn));
for i = 1:numel(Info.Stat)
    for j = setdiff(1:numel(Info.Stat), i)
        text(j, i, sprintf('%2.0f%%', Comp.(Var).BayesStat.normevidence(i,j)*100));
        plot([0,numel(Info.Stat)+1], repmat(i-0.5,1,2), 'k-', 'LineWidth', lw)
        plot(repmat(j-0.5,1,2), [0,numel(Info.Stat)+1], 'k-', 'LineWidth', lw)
    end
end
set(gca, 'XTick', 1:numel(Info.Stat), 'XTickLabel', Info.Stat); 
set(gca, 'YTick', 1:numel(Info.Stat), 'YTickLabel', Info.Stat, 'XTickLabelRotation', 45);
set(gca, 'Box', 'On', 'LineWidth', lw); axis('square');
title({'Normalized evidence ratios','',''});

%
subplot(2,2,3);
imagesc(Comp.(Var).BayesIntg.normevidence); hold('on'); colormap(flipud(autumn));
for i = 1:numel(Info.Intg)
    for j = setdiff(1:numel(Info.Intg), i)
        text(j, i, sprintf('%2.0f%%', Comp.(Var).BayesIntg.normevidence(i,j)*100));
        plot([0,numel(Info.Intg)+1], repmat(i-0.5,1,2), 'k-', 'LineWidth', lw)
        plot(repmat(j-0.5,1,2), [0,numel(Info.Intg)+1], 'k-', 'LineWidth', lw)
    end
end
set(gca, 'XTick', 1:numel(Info.Intg), 'XTickLabel', Info.Intg);
set(gca, 'YTick', 1:numel(Info.Intg), 'YTickLabel', Info.Intg, 'XTickLabelRotation', 45);
set(gca, 'Box', 'On', 'LineWidth', lw); axis('square');

%% PLOT THE QUALITY OF FIT
%  =======================

% Define Bayesian Information Criterion function
BIC = @(MSE,k,n) (n .* log(MSE)) + (k .* log(n));

% Get the number of observations
pattofit = Info.FitPatLen(nPat);
n = sum(sum(sum(~isnan(Data(pattofit,:,:)))));

% Make a BIC map based on k and MSE
kgrid = -1:1:11;
if     strcmp(d, 'SquiresScience1976'), LLHgrid = linspace(0.4, 1.6, 1000);
elseif strcmp(d, 'KolossaFIHN2013'),    LLHgrid = linspace(0.1, 1.2, 1000);
end
BICmat = NaN(numel(LLHgrid), numel(kgrid));
for i = 1:numel(kgrid)
    BICmat(:,i) = BIC(LLHgrid, kgrid(i), n);
end

% Prepare stuffs
dotsmrk = {'^', 'v', 'd', 'o'};
dotscol = {'k', 'k', 'k', 'w', 'k'};
lgd = [];
figure('Units', 'Normalized', 'Position', [0, 0, 0.25, 0.45]);

% Plot the BIC grid
imagesc(kgrid, LLHgrid, BICmat); hold('on');
%surf(kgrid, LLHgrid, BICmat, 'EdgeColor', 'None', 'FaceAlpha', 0.5); hold('on');
contour(kgrid, LLHgrid, BICmat, 'EdgeColor', 'k', 'ShowText', 'On');

% For each Bayesian observer
for s = 1:nStat
    for i = nIntg:-1:1
        
        % Plot the location of that model on the BIC space
        lgd(s,:) = plot(Best.k(i,s), Best.MSE(i,s), dotsmrk{s}, ...
            'MarkerFaceColor', dotscol{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
        text(Best.k(i,s)+0.3, Best.MSE(i,s), sprintf('Learning %s with a %s integration', ...
            Info.Stat{s}, Info.Intg{i}), 'HorizontalAlignment', 'Left', 'FontSize', fs/2);
%         lgd(s,:) = plot3(Best.k(i,s), Best.MSE(i,s), Best.BIC(i,s), dotsmrk{s}, ...
%             'MarkerFaceColor', dotscol{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
%         text(Best.k(i,s)+0.3, Best.MSE(i,s), Best.BIC(i,s), sprintf('Learning %s with a %s integration', ...
%             Info.Stat{s}, Info.Intg{i}), 'HorizontalAlignment', 'Left', 'FontSize', fs/2);
    end
end

% For each descriptive model
for m = 1:numel(Info.Desc)
    
    % Plot the location of that model on the BIC space
    lgd(s,:) = plot(Desc.k(m), Desc.MSE(m), dotsmrk{end}, ...
        'MarkerFaceColor', dotscol{end}, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    text(Desc.k(m)+0.3, Desc.MSE(m), sprintf('%s model', Info.Desc{m}), ...
        'HorizontalAlignment', 'Left', 'FontSize', fs/2);
end

% Customize the axes
cbr = colorbar('SouthOutside'); cbr.LineWidth = lw; cbr.Label.String = 'BIC'; cbr.Label.FontSize = fs;
colormap(parula); axis('xy'); axis([min(kgrid)+0.5, max(kgrid)-0.5, min(LLHgrid), max(LLHgrid)]);
set(gca, 'LineWidth', lw, 'FontSize', fs, 'Box', 'On', 'Layer', 'Top');
set(gca, 'XTick', min(kgrid):1:max(kgrid), 'YTick', min(LLHgrid):0.1:max(LLHgrid));
legend(lgd, Info.Stat, 'Location', 'NorthEastOutside');
xlabel('Number of free parameters (k)', 'FontSize', fs);
ylabel('Mean squared error (MSE)', 'FontSize', fs);
title({dataset{1}, sprintf('Value of BIC given n = %i', n), ''}, 'FontSize', fs);

%% PLOT THE REGRESSION BETWEEN SURPRISE AND REAL DATA
% ===================================================

% Make the data x predictions plots
figure('Units', 'Normalized', 'Position', [0.28 0 0.4 0.45]);
x = min(Data(:)):0.001:max(Data(:));

% For each model
for i = 1:nIntg
    for s = 1:nStat
        subplot(nStat,nIntg,i+(nIntg*(s-1))); cla;
        
        % Plot the regression
        plot(x, x, '-', 'color', g, 'LineWidth', lw*2); hold('on');
        
        % For each global probability block
        for p = 1:nP
            
            % Plot the dots corresponding to each pattern
            plot(Data(pattofit,:,p), Best.Simu{i,s}(pattofit,:,p), 'o', 'MarkerSize', ...
                6, 'MarkerFaceColor', col(p,:), 'MarkerEdgeColor', 'k', 'LineWidth', 0.25); hold('on');
        end
        
        % Customize the axes
        axis('square'); xlim([-1, 6]); ylim(xlim);
        set(gca, 'XTick', -1:1:6, 'XGrid', 'On', 'YTick', -1:1:6, 'YGrid', 'On');
        set(gca, 'LineWidth', lw, 'FontSize', fs, 'Box', 'On', 'Layer', 'Bottom', 'TickDir', 'In');
        if     i == 1,     ylabel({['Learning ', lower(Info.Stat{s})], '', 'Fitted surprise'}, 'FontWeight', 'Normal', 'FontSize', fs); end
        if     s == 1,     title({[Info.Intg{i}, ' integration'], ''}, 'FontWeight', 'Normal', 'FontSize', fs);
        elseif s == nStat, xlabel('P300 amplitude', 'FontWeight', 'Normal', 'FontSize', fs); end
    end
end

%% PLOT THE TREES
%  ==============

if     strcmpi(d, 'SquiresScience1976'), axlim = [3.5, 11, -0.5, 5.5];
elseif strcmpi(d, 'KolossaFIHN2013'),    axlim = [3, 9, -0.5, 5.5]; end

% Plot the "data trees" (i.e. the one we are trying to fit)
figure('Units', 'Normalized', 'Position', [0.7, 0.52, 0.35, 0.55]); lgd = [];
for p = 1:nP
    treepos = (nPat/2)*(p-1)+2;
	lgd(p,:) = PlotTree2(Data(:,:,p), nPat, [1 1 1; 0 0 0], lw/2, treepos, 6, 0.08, col(p,:), lw);
end
set(gca, 'XTick', [], 'YTick', -0.5:0.5:5.5, 'YGrid', 'On', 'YMinorGrid', 'Off', 'Layer', 'Bottom');
set(gca, 'LineWidth', lw, 'Box', 'Off');
axis(axlim);
legend(lgd, cellfun(@(x) sprintf('p(o) = %1.2f', x), num2cell(pX), ...
    'UniformOutput', false), 'Location', 'SouthEast', 'FontSize', fs);
ylabel({'P300 averaged amplitude', ''}, 'FontSize', fs);
title(dataset, 'FontWeight', 'Bold', 'FontSize', fs);

% Plot the theoretical trees (those with the best parameter fit)
iter = 1;
for s = [3,1,2]
    for i = [1,2]
        figure('Units', 'Normalized', 'Position', [0.7, 0.52, 0.35, 0.55]); iter = iter + 1;
        for p = 1:nP
            treepos = (nPat/2)*(p-1)+2;
            lgd(p,:) = PlotTree2(Best.Simu{i,s}(:,:,p), nPat, [1 1 1; 0 0 0], lw/2, treepos, 6, 0.08, col(p,:), lw);
        end
        set(gca, 'XTick', [], 'YTick', -0.5:0.5:5.5, 'YGrid', 'On', 'YMinorGrid', 'Off', 'Layer', 'Bottom');
        set(gca, 'LineWidth', lw, 'Box', 'Off');
        axis(axlim);
        ylabel({'Fitted surprise', ''}, 'FontSize', fs);
        title({sprintf('Learning %s with a %s memory', lower(Info.Stat{s}), ...
            lower(Info.Intg{i})), ''}, 'FontWeight', 'Bold', 'FontSize', fs);
    end
end
