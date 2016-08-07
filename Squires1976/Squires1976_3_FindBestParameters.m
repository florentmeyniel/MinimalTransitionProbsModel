% Return best parameters of the fit and BIC (AIC) values.
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

% Which quantity to fit?
q = 1; % surprise
%q = 2; % update

% Which fitting procedure to use?
w = 1; % not weighted
%w = 2; % weighed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the functions
addpath(genpath('Squires1976'));
try cd('Squires1976'); catch, end

%% LOAD DATA
%  =========

% Load the results
load(sprintf('Squires1976_FitModels_%s.mat', d));

% Add information
Info.Quty = Info.Quty{q};
Info.IsWeighted = logical(w-1);
if ~Info.IsWeighted
    Info = rmfield(Info, 'WeightFun');
end

% Load data
[Data, nPat] = Squires1976_PrepareData(sprintf('Data%s.mat', d));

%% FIND THE BEST MODELS
%  ====================

% Get the models to plot
modlist = Info.Models;
modtoplot = modlist(~strcmpi(Info.Models(:,1), 'Descriptive'),:);
[~,idx] = unique(modtoplot(:,1)); % avoid sorting
stat = modtoplot(sort(idx),1);
nStat = numel(stat);
[~,idx] = unique(modtoplot(:,2)); % avoid sorting
intg = modtoplot(sort(idx),2);
nIntg = numel(intg);

% Define Bayesian Information Criterion function
BIC = @(MSE,k,n) (n .* log(MSE)) + (k .* log(n));

% Define the Akaike Information Criterion function
AIC = @(RSS,k,n) 2*k + n*log(RSS);

% Get the number of observations
pattofit = Info.FitPatLen(nPat);
n = sum(sum(sum(~isnan(Data(pattofit,:,:)))));

% Prepare outputs
nMod = numel(Fit.MSE);
Best.MSE   = NaN(nIntg+1,nStat);
Best.Idx   = NaN(nIntg+1,nStat);
Best.k     = NaN(nIntg+1,nStat);
Best.BIC   = NaN(nIntg+1,nStat);
Best.AIC   = NaN(nIntg+1,nStat);
Best.R2    = NaN(nIntg+1,nStat);
Best.Param = NaN(nIntg+1,nStat);
Best.Grid  = cell(nIntg+1,nStat);
Best.Simu  = cell(nIntg+1,nStat);

% For each statistic
for s = 1:nStat

    % For each integration type
    for i = 0:nIntg % i = 0 => perfect integration

        % Find the model's specifications
        if i == 0
            m = find(strcmpi(Info.Models(:,1), stat{s}) + ... % statistic
                     strcmpi(Info.Models(:,2), 'Windowed') == 2); % integration
        else
            m = find(strcmpi(Info.Models(:,1), stat{s}) + ... % statistic
                     strcmpi(Info.Models(:,2), intg{i}) == 2); % integration
        end

        % Find the minimum error
        if i == 0
            Best.MSE(1,s) = Fit.MSE{m}(end,q,w);
            Best.Idx(1,s) = numel(Fit.MSE{m}(:,q,w));            
        else
            [Best.MSE(i+1,s), Best.Idx(i+1,s)] = min(Fit.MSE{m}(:,q,w));
        end
        
        % Compute the BIC
        if i == 0
            Best.k(1,s) = 0; % they do not have any free parameter
        else
            Best.k(i+1,s) = 1; % they all have one single free parameter
        end
        Best.k(i+1,s) = Best.k(i+1,s) + 2; % offset and scaling
        Best.BIC(i+1,s) = BIC(Best.MSE(i+1,s), Best.k(i+1,s), n);
        Best.AIC(i+1,s) = AIC(Best.MSE(i+1,s), Best.k(i+1,s), n);
        
        % Find the corresponding R2
        Best.R2(i+1,s) = Fit.R2{m}(Best.Idx(i+1,s),q,w);

        % Find the associated parameter
        if i == 0
            Best.Grid{1,s} = Spec.ParamGrid{m}(Best.Idx(1,s));
        else
            Best.Grid{i+1,s} = Spec.ParamGrid{m};
        end
        Best.Param(i+1,s) = Spec.ParamGrid{m}(Best.Idx(i+1,s));

        % Get the corresponding model's predictions
        Best.Simu{i+1,s} = Fit.ScSurp{m}(:,:,:,Best.Idx(i+1,s),q,w);
    end
end

% FOR THE DESCRIPTIVE MODEL
% ~~~~~~~~~~~~~~~~~~~~~~~~~

% Prepare outputs
Desc.MSE   = NaN(1,2);
Desc.Idx   = NaN(1,2);
Desc.k     = NaN(1,2);
Desc.BIC   = NaN(1,2);
Desc.AIC   = NaN(1,2);
Desc.R2    = NaN(1,2);
Desc.Grid  = NaN(1,2);
Desc.Param = NaN(1,2);
Desc.Simu  = cell(1,2);

% Only the surprise is available
if q == 1

    % For both Squires and Kolossa models
    for m = 1:2

        % Find the model's specifications
        mm = size(modtoplot,1)+m;

        % Get the minimum error
        Desc.MSE(m) = Fit.MSE{mm}(1,q,w);
        Desc.Idx(m) = 1;

        % Compute the BIC
        if     strfind(Info.Models{mm,2}, 'Squires'), Desc.k(m) = 3;
        elseif strfind(Info.Models{mm,2}, 'Kolossa'), Desc.k(m) = 6; end
        Desc.k(m) = Desc.k(m) + 2; % offset and scaling
        Desc.BIC(m) = BIC(Desc.MSE(m), Desc.k(m), n);
        Desc.AIC(m) = AIC(Desc.MSE(m), Desc.k(m), n);

        % Get the R2
        Desc.R2(m) = Fit.R2{mm}(1,q,w);

        % Parameters not specified
        Desc.Grid(m)  = NaN;
        Desc.Param(m) = NaN;

        % Get the model's predictions
        Desc.Simu{m} = Fit.ScSurp{mm}(:,:,:,1,q,w);
    end
end

%% PERFORM MODELS COMPARISON
%  =========================

% Deltas (difference from the best model)
BICs = [Best.BIC(:); Desc.BIC(:)];
Comp.BIC.delta = BICs - min(BICs);
AICs = [Best.AIC(:); Desc.AIC(:)];
Comp.AIC.delta = AICs - min(AICs);

% Pairwise difference and ratios between models
eve = repmat(BICs, [1,numel(BICs)]);
Comp.BIC.BF.div = eve ./ eve';
Comp.BIC.BF.sub = eve -  eve';
eve = repmat(AICs, [1,numel(AICs)]);
Comp.AIC.BF.div = eve ./ eve';
Comp.AIC.BF.sub = eve -  eve';

% For each measured variable
Var = fieldnames(Comp);
for v = 1:numel(Var)
    
    % Models' weights
	Comp.(Var{v}).All.weights   = exp((-1/2) .* Comp.(Var{v}).delta) ./ sum(exp((-1/2) .* Comp.(Var{v}).delta));
    Comp.(Var{v}).BayesAll.weights = exp((-1/2) .* Comp.(Var{v}).delta(1:end-2)) ./ sum(exp((-1/2) .* Comp.(Var{v}).delta(1:end-2)));

    % Reshape such that the matrix as the same format as the previous variables
    Comp.(Var{v}).BayesAll.weights = reshape(Comp.(Var{v}).BayesAll.weights, size(Best.BIC));
    
    % Weights of models' families
    Comp.(Var{v}).BayesStat.weights = sum(Comp.(Var{v}).BayesAll.weights,1);
    Comp.(Var{v}).BayesIntg.weights = sum(Comp.(Var{v}).BayesAll.weights,2);
    
    % For each group of models
    Cat = fieldnames(Comp.(Var{v}));
    Cat = Cat(3:end);
    for c = 1:numel(Cat)
        
        % Compare weights
        eve = repmat(Comp.(Var{v}).(Cat{c}).weights(:), [1, numel(Comp.(Var{v}).(Cat{c}).weights)]);
        Comp.(Var{v}).(Cat{c}).evidence     = eve ./ eve';
        Comp.(Var{v}).(Cat{c}).normevidence = eve ./ (eve + eve');
    end
end

%% SAVE THE RESULTS
%  ================

% Update the info fields
Info.Stat = stat;
Info.Intg = vertcat('Perfect', intg);
Info.Desc = Info.Models(end-1:end,2);
Info.Parameters = {'Window size', 'Leak strength', 'p(change)'};
Spec.Models = Info.Models;
Info = rmfield(Info, 'Models');

% Save the file
save(sprintf('Squires1976_BestParameters_%s.mat', d), '-v7.3', 'Info', 'Spec', 'Best', 'Desc', 'Comp');

%% PLOT THE QUALITY OF FIT OVER PARAMETERS
%  =======================================

% Font size and line width for the plots
fs = 12;
lw = 1.5;
intgcol = lines(nIntg);

% Options for the inset plot
inset_w = 0.305; % fraction of subplot width
inset_h = 0.215; % fraction of subplot width
inset_gapw  = [0.00, 0.00, 0.00]; % fraction of subplot width
inset_gaph  = [0.00, 0.00, 0.00]+0.035; % fraction of subplot width;
inset_rectw = [0.02, 0.03, 0.005]; % fraction of the x-axis length
inset_recth = [0.07, 0.09, 0.04]; % fraction of the y-axis length

% Figure options
figure('Units', 'Normalized', 'Position', [0.15 0.35 0.75 0.25]);
xt = {[1 2 4 8 15 30 55 100 200], ...
      [1 2 4 8 15 30 55 100 200], ...
      [0.001 0.003 0.01 0.03 0.1 0.3 1]};

% Plot each integration type
for i = 1:nIntg
    subplot(1,nIntg,i); lgd = [];

    % For each statistic
    for s = 1:nStat

        % Find the model
        m = find(strcmpi(Spec.Models(:,1), stat{s}) + ... % statistic
                 strcmpi(Spec.Models(:,2), intg{i}) == 2); % integration

        % Plot the error
        lgd(s,:) = plot(Spec.ParamGrid{m}, Fit.MSE{m}(:,q,w), '.-', 'Color', intgcol(s,:), ...
            'MarkerSize', lw*10, 'LineWidth', lw, 'MarkerFaceColor', 'None'); hold('on');

        % Plot the best model
        plot(Best.Param(i+1,s), Best.MSE(i+1,s), 'o', 'Color', intgcol(s,:), ...
            'MarkerSize', lw*7,  'LineWidth', lw);
        text(Best.Param(i+1,s), Best.MSE(i+1,s), sprintf('%g', Best.Param(i+1,s)), ...
            'VerticalAlignment', 'Bottom', 'Color', intgcol(s,:), 'FontWeight', 'Bold', 'FontSize', fs);
    end
    
    % Axes option
    ylimmax = max(cellfun(@(x) max(x(:)), Fit.MSE));
    axis([Spec.ParamGrid{m}(1), Spec.ParamGrid{m}(end), 0, ylimmax]);
    set(gca, 'Box', 'Off', 'Layer', 'Bottom', 'LineWidth', lw, 'FontSize', fs);
    set(gca, 'XScale', 'log', 'XTick', xt{i}, 'XGrid', 'On', 'XMinorTick', 'Off', 'XMinorGrid', 'Off');
    set(gca, 'XTickLabel', get(gca, 'XTick'));
    set(gca, 'YTick', linspace(0, ylimmax, 10));
    set(gca, 'YTickLabel', num2str(get(gca,'YTick')', '%1.2f'));

    % Labels
    if i == 1, ylabel('Mean squared error', 'FontSize', fs);
    elseif i == 3, l = legend(lgd, stat, 'Location', 'SouthEast', 'Box', 'On'); end
    title({[intg{i}, ' integration'], ''}, 'FontWeight', 'Bold', 'FontSize', fs);
    xlabel(Info.Parameters{i}, 'FontSize', fs);

    % Get the inset dimensions
    inset_xorigin = mean(Best.Param(i+1,[1,3]));
    inset_width   = inset_rectw(i) * (Spec.ParamGrid{m}(end) - Spec.ParamGrid{m}(1) + 1);
    inset_xlim    = [inset_xorigin - inset_width, inset_xorigin + inset_width];
    inset_yorigin = mean(Best.MSE(i+1,[1,3]));
    inset_heigth  = inset_recth(i) * diff(ylim);
    inset_ylim    = [inset_yorigin - inset_heigth, inset_yorigin + inset_heigth];
    if inset_xlim(1) < 0, inset_xlim(1) = 0; end

%     % Plot a rectangle
%     fill([inset_xlim(1), inset_xlim(2), inset_xlim(2), inset_xlim(1)], ...
%          [inset_ylim(1), inset_ylim(1), inset_ylim(2), inset_ylim(2)], ...
%           '--', 'EdgeColor', repmat(128/255,1,3), 'FaceColor', 'None', 'LineWidth', lw/2);

    % Create axes for the inset plot
    pos = get(gca, 'Position');
    pos = [pos(1) + pos(3) - (pos(3) * inset_w) - (pos(3) * inset_gapw(i)), ...
           pos(2) + pos(4) - (pos(4) * inset_h) - (pos(4) * inset_gaph(i)), ...
           (pos(3) * inset_w), (pos(4) * inset_h)];
    ax = axes('Position', pos, 'Box', 'On', 'LineWidth', lw);

    % Plot the data in the inset plot
    for p = 1:numel(Spec.ParamGrid{m})
        plot(Spec.ParamGrid{m}(p), Fit.MSE{m}(p,q,w), 'k:'); hold('on');
    end

    % For each statistic
    for s = 1:nStat
        
        % Find the model
        m = find(strcmpi(Spec.Models(:,1), stat{s}) + ... % statistic
                 strcmpi(Spec.Models(:,2), intg{i}) == 2); % integration
             
        % Plot the error
        plot(squeeze(Spec.ParamGrid{m}), Fit.MSE{m}(:,q,w), '.-', 'Color', intgcol(s,:), ...
            'MarkerSize', lw*10, 'LineWidth', lw, 'MarkerFaceColor', 'None'); hold('on');
        
        % Plot the best model
        plot(Best.Param(i+1,s), Best.MSE(i+1,s), 'o', 'Color', intgcol(s,:), ...
            'MarkerSize', lw*7,  'LineWidth', lw);
    end

    % Scale the inset plot such as it acts as a zoom in
    axis([inset_xlim, inset_ylim]);
    set(gca, 'Box', 'On', 'Layer', 'Bottom', 'LineWidth', lw/2, 'FontSize', fs*0.5);
    set(gca, 'XTick', Spec.ParamGrid{m}, 'XTickLabelRotation', 45, 'XGrid', 'On', 'XAxisLocation', 'Top');
    set(gca, 'YTick', linspace(inset_ylim(1), inset_ylim(2), 5), 'YGrid', 'On', 'YAxisLocation', 'Right');
    set(gca, 'YTickLabel', num2str(get(gca,'YTick')', '%1.3f'));
end
