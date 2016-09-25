% Answer to one of the reviewers:
%   - identifiability analysis
%   - cross-validation analysis

% Copyright 2016 Florent Meyniel & Maxime Maheu

%% INITIALIZATION
%  ==============

% Clear the place
clear; close('all'); clc;

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

% Load data
[data, nPat, p1] = Squires1976_PrepareData(sprintf('Data%s.mat', d));
nP = numel(p1);

% Reorder the blocks
[~,I] = sort(p1, 2, 'descend');
data = data(:,:,I);
p1 = p1(I);

% Keep only the data from the longest patterns
data = squeeze(data(end,:,:));
data = data(:);

% Load simulations corresponding to all the parameters
A = load(sprintf('Squires1976_FitModels_%s.mat', d));

% Load simulations corresponding to the best parameter
B = load(sprintf('Squires1976_BestParameters_%s.mat', d));

%% PERFORM THE CROSS-VALIDATION
%  ============================

% Note that for models with a perfect integration, we pick the model with a
% windowed integration and a window parameter at 200 (because longest
% sequences used for simulations are composed of 200 stimuli).

% Define models
stat = B.Info.Stat;
intg = B.Info.Intg;
tmp = repmat(intg, [numel(stat),1]);
tmp = tmp([1:4:end,2:4:end,3:4:end,4:4:end]);
models = [repmat(stat, [numel(intg),1]), tmp];
nmods = size(models,1);

% Prepare outputs
ID_epsilon = NaN(nmods,nmods,numel(data));
ID_bestpar = NaN(nmods,nmods,numel(data));
CV_epsilon = NaN(nmods,numel(data));
CV_bestpar = NaN(nmods,numel(data));

% Loop over models
for m = 1:nmods

    % Get the model
    fprintf('* Model %i/%i: "%s" statistic with a "%s" integration...\n', ...
        m, nmods, models{m,1}, models{m,2});

    % Get the indices corresponding to the current model
    s = find(strcmpi(models{m,1}, B.Info.Stat));
    if   strcmpi(models{m,2}, 'Perfect'), i = find(strcmpi('Windowed',  B.Info.Intg));
    else                                  i = find(strcmpi(models{m,2}, B.Info.Intg)); end
    midx = intersect(find(strcmpi(A.Info.Models(:,1), stat{s})), find(strcmpi(A.Info.Models(:,2), intg{i})));

    % Get the simulation (i.e. trees) corresponding to the best parameter
    if strcmpi(models{m,2}, 'Perfect')
        simu = B.Best.Simu{1,s};
        bidx = numel(A.Spec.ParamGrid{midx});
    else
        simu = B.Best.Simu{i,s};
        bidx = B.Best.Idx(i,s);
    end

    % Back to the surprise scale (the trees was saved in the scaled format)
    beta1 = A.Fit.Beta{midx}(1,bidx,q,w);
    beta0 = A.Fit.Beta{midx}(2,bidx,q,w);
    simu = (simu - beta0) ./ beta1;

    % Keep only the longest patterns
    simu = squeeze(simu(end,:,:));
    n = numel(simu);
    simu = reshape(simu, [n,1]);

    % Perform a leave-one-out procedure
    for j = 1:n
        fprintf('  - leave-one-out %i/%i,\n', j, n);

        % Remove one observation
        idx = setdiff(1:n,j);

        % Loop (again) over models
        for mm = 1:nmods

            % Get the indices corresponding to the current model
            ss = find(strcmpi(models{mm,1}, B.Info.Stat));
            if strcmpi(models{mm,2}, 'Perfect'), ii = find(strcmpi('Windowed',   B.Info.Intg));
            else                                 ii = find(strcmpi(models{mm,2}, B.Info.Intg)); end
            mmidx = intersect(find(strcmpi(A.Info.Models(:,1), stat{ss})), find(strcmpi(A.Info.Models(:,2), intg{ii})));

            % Get parameter's grid
            if strcmpi(models{mm,2}, 'Perfect'), params = A.Spec.ParamGrid{mmidx}(end);
            else                                 params = A.Spec.ParamGrid{mmidx}; end
            np = numel(params);

            % Loop over parameters
            ID_MSE  = NaN(1,np);
            ID_Beta = NaN(2,np);
            ID_Surp = NaN(n,np);
            if mm == m
                CV_MSE  = NaN(1,np);
                CV_Beta = NaN(2,np);
                CV_Surp = NaN(n,np);
            end
            for p = 1:np

                % For perfect models, find the last (perfect) integration parameter
                if np == 1
                    p = numel(A.Spec.ParamGrid{1});
                    pi = 1;
                else pi = p;
                end

                % Get the corresponding (scaled) tree
                pred = A.Fit.ScSurp{mmidx}(:,:,:,p,q,w);

                % Back to the surprise scale (it was saved in the scaled format)
                beta1 = A.Fit.Beta{mmidx}(1,p,q,w);
                beta0 = A.Fit.Beta{mmidx}(2,p,q,w);
                pred = (pred - beta0) ./ beta1;

                % Keep only the longest patterns
                pred = squeeze(pred(end,:,:));
                pred = reshape(pred, [n,1]);

                % Save the predictions (in -log(p) scale for latter)
                ID_Surp(:,pi) = pred;

                % Fit the reduced dataset (i.e. the predictions of model
                % "m" minus 1 observation) with the model "mm" (with a
                % parameter "p")
                warning('off');
                [ID_MSE(pi),~,~,~,~,~,ID_Beta(:,pi)] = Squires1976_FitModel(pred(idx), simu(idx));
                warning('on');

                % Fit the reduced dataset (i.e. the Squires' data minus
                % one observation) with the model "m" (with a parameter "p")
                if mm == m
                    CV_Surp(:,pi) = pred;
                    warning('off');
                    [CV_MSE(pi),~,~,~,~,~,CV_Beta(:,pi)] = Squires1976_FitModel(pred(idx), data(idx));
                    warning('on');
                end
            end

            % IDENTIFIABILITY
            % ~~~~~~~~~~~~~~~
            % A particular model (model mm) is fitted on the predictions of an
            % other model (model m, with the best fitted parameter from the fit
            % of Squires' data). The prediction of this model (model mm, with
            % the best fitted parameter) for the leave one out observation is
            % compared to the prediction of the simulated model (model m).

            % Get the best parameter
            [~,bpi] = min(ID_MSE);
            if     numel(ID_MSE) == 1, ID_bestpar(m,mm,j) = A.Spec.ParamGrid{mmidx}(end);
            elseif numel(ID_MSE) >  1, ID_bestpar(m,mm,j) = A.Spec.ParamGrid{mmidx}(bpi); end

            % Get the corresponding prediction for that particular
            % leave out observation
            pred = (ID_Surp(:,bpi) .* ID_Beta(1,bpi)) + ID_Beta(2,bpi);

            % Get the prediction for the observation that was left out and
            % measure its distance to the data
            ID_epsilon(m,mm,j) = simu(j) - pred(j);
        end

        % CROSS VALIDATION
        % ~~~~~~~~~~~~~~~~
        % A particular model is fitted on Squires' data minus one
        % observation. The prediction of this model (with the best fitted
        % parameter) for the leave out observation is compared to the real
        % observation from Squires et al.

        % Get the best parameter
        [~,bpi] = min(CV_MSE);
        if     numel(CV_MSE) == 1, CV_bestpar(m,j) = A.Spec.ParamGrid{midx}(end);
        elseif numel(CV_MSE) >  1, CV_bestpar(m,j) = A.Spec.ParamGrid{midx}(bpi); end

        % Get the corresponding prediction for that particular
        % leave out observation
        pred = (CV_Surp(:,bpi) .* CV_Beta(1,bpi)) + CV_Beta(2,bpi);

        % Get the prediction for the observation that was left out and
        % measure its distance to the data
        CV_epsilon(m,j) = data(j) - pred(j);
    end
    fprintf('\n');
end

%% PLOT THE RESULTS OF THE IDENTIFIABILITY ANALYSIS
%  ================================================

% Indices of models to look at
perfidx = 1:3;
leakidx = 7:9;

% Average of error in specific models
Eavg_AtoB{1} = mean(abs(ID_epsilon(leakidx,leakidx,:)),3); % leak to leak
Eavg_AtoB{2} = mean(abs(ID_epsilon(leakidx,perfidx,:)),3); % leak to perfect
Eavg_AtoB{3} = mean(abs(ID_epsilon(perfidx,leakidx,:)),3); % perfect to leak
Eavg_AtoB{4} = mean(abs(ID_epsilon(perfidx,perfidx,:)),3); % perfect to perfect

% Standard deviation of error in specific models
Estd_AtoB{1} = std(abs(ID_epsilon(leakidx,leakidx,:)),[],3); % leak to leak
Estd_AtoB{2} = std(abs(ID_epsilon(leakidx,perfidx,:)),[],3); % leak to perfect
Estd_AtoB{3} = std(abs(ID_epsilon(perfidx,leakidx,:)),[],3); % perfect to leak
Estd_AtoB{4} = std(abs(ID_epsilon(perfidx,perfidx,:)),[],3); % perfect to perfect

% Useful variables for the plot
from = {'leaky', 'leaky', 'perfect', 'perfect'};
to   = {'leaky', 'perfect', 'leaky', 'perfect'};
pos  = [-2/9, 0, 2/9];
cols = lines(3);

% Loop over pairs of integration types
figure('Name', 'Identifiability', 'Position', [0.1859 0.1683 0.6286 0.6250]);
for i = 1:4
    subplot(2,2,i);

    % Plot the bars
    bar(1:numel(stat), Eavg_AtoB{i}, 1, 'LineWidth', 1); hold('on');

    % Plot the error bars and the text labels
    for m = 1:3
        for mm = 1:3
            M = Eavg_AtoB{i}(m,mm);
            S = Estd_AtoB{i}(m,mm) ./ sqrt(n);
            plot(repmat(m + pos(mm),1,2), [M-S, M+S], 'k-', 'LineWidth', 1);
            text(m + pos(mm), 0.03, sprintf('%1.2f', M), 'Color', 'k');
        end
    end

    % Customize the plot
    ylim([0,0.5]);
    set(gca, 'XTick', 1:numel(stat), 'XTickLabel', stat, 'XTickLabelRotation', 0);
    set(gca, 'Box', 'Off', 'LineWidth', 1, 'FontSize', 14);
    lgd = legend(stat, 'Box', 'Off', 'Location', 'North', 'Orientation', 'Horizontal');
    colormap(cols);

    % Add labels
    xlabel(sprintf('From %s ...', from{i}));
    title(sprintf('... to %s ...', to{i}), 'FontWeight', 'Normal');
    ylabel({'Absolute distance averaged', 'over simulations (+/- SEM)'});
end

%% PLOT THE RESULTS OF THE CROSS-VALIDATION ANALYSIS
%  =================================================

% Average of error in specific models
Eavg_AtoB{1} = mean(abs(CV_epsilon(perfidx,:)),2); % leak
Eavg_AtoB{2} = mean(abs(CV_epsilon(leakidx,:)),2); % perfect

% Standard deviation of error in specific models
Estd_AtoB{1} = std(abs(CV_epsilon(perfidx,:)),[],2); % leak
Estd_AtoB{2} = std(abs(CV_epsilon(leakidx,:)),[],2); % perfect

% Average of best parameters in specific models
Pavg_AtoB{1} = mean(CV_bestpar(perfidx,:),2); % leak
Pavg_AtoB{2} = mean(CV_bestpar(leakidx,:),2); % perfect

% Standard deviation of best parameters in specific models
Pstd_AtoB{1} = std(CV_bestpar(perfidx,:),[],2); % leak
Pstd_AtoB{2} = std(CV_bestpar(leakidx,:),[],2); % perfect

% Useful variables for the plot
intgt = intg([1,3]);

% Loop over models
figure('Name', 'Cross-validation', 'Position', [0.2698 0.2508 0.4609 0.4600]);
for i = 1:2
    subplot(1,2,i);
    for m = 1:3

        % Plot the bars
        M = Eavg_AtoB{i}(m);
        bar(m, abs(M), 'FaceColor', cols(m,:), 'LineWidth', 1); hold('on');

        % Plot the error bars
        S = Estd_AtoB{i}(m) ./ sqrt(n);
        plot(repmat(m,1,2), [M-S, M+S], 'k-', 'LineWidth', 1);

        % Plot the best fitted parameters
        text(m, 0.1, sprintf('%1.0f', Pavg_AtoB{i}(m)), 'Color', 'k');
        text(m, 0.05, sprintf('(%.2f)', Pstd_AtoB{i}(m)), 'Color', 'k', 'FontAngle', 'Italic');
    end

    % Customize the plot
    axis([0,4,0,1.2]);
    set(gca, 'Box', 'Off', 'LineWidth', 1, 'FontSize', 14);

    % Add labels
    set(gca, 'XTick', 1:3, 'XTickLabel', stat);
    title(sprintf('Models with a %s integration', lower(intgt{i})));
    ylabel('Absolute distance averaged over simulations (+/- SEM)');
end
