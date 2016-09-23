% Answer to one of the reviewers: cross-validation analysis

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
[Data, nPat, p1] = Squires1976_PrepareData(sprintf('Data%s.mat', d));
nP = numel(p1);

% Reorder the blocks
[~,I] = sort(p1, 2, 'descend');
Data = Data(:,:,I);
p1 = p1(I);

% Keep only the data from the longest patterns
Data = squeeze(Data(end,:,:));
Data = Data(:);

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
epsilon = NaN(nmods,nmods,numel(Data));
bestp = NaN(nmods,nmods,numel(Data));
R2within = NaN(nmods,nmods,numel(Data));
R2generalized = NaN(nmods,nmods,numel(Data));

% Loop over models
for m = 1:nmods
    
    % Get the model
    fprintf('* Model %i/%i: "%s" statistic with a "%s" integration...\n', ...
        m, nmods, models{m,1}, models{m,2});

    % Get the corresponding indices
    s = find(strcmpi(models{m,1}, B.Info.Stat));
    if   strcmpi(models{m,2}, 'Perfect'), i = find(strcmpi('Windowed',  B.Info.Intg));
    else                                  i = find(strcmpi(models{m,2}, B.Info.Intg)); end
    midx = intersect(find(strcmpi(A.Info.Models(:,1), stat{s})), find(strcmpi(A.Info.Models(:,2), intg{i})));

    % Get the corresponding best simulation
    if strcmpi(models{m,2}, 'Perfect')
        simu = B.Best.Simu{1,s};
        bidx = numel(A.Spec.ParamGrid{midx});
    else
        simu = B.Best.Simu{i,s};
        bidx = B.Best.Idx(i,s);
    end
    
    % Back to the surprise scale (it was saved in the scaled format)
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
        loo = simu(idx);
    
        % Loop (again) over models
        for mm = 1:nmods
            
            % Get the corresponding indices
            ss = find(strcmpi(models{mm,1}, B.Info.Stat));
            if strcmpi(models{mm,2}, 'Perfect'), ii = find(strcmpi('Windowed',   B.Info.Intg));
            else                                 ii = find(strcmpi(models{mm,2}, B.Info.Intg)); end
            mmidx = intersect(find(strcmpi(A.Info.Models(:,1), stat{ss})), find(strcmpi(A.Info.Models(:,2), intg{ii})));
            
            % Get parameters
            if strcmpi(models{mm,2}, 'Perfect'), params = A.Spec.ParamGrid{mmidx}(end);
            else                                 params = A.Spec.ParamGrid{mmidx}; end
            np = numel(params);

            % Loop over parameters
            R2   = NaN(1,np);
            MSE  = NaN(1,np);
            beta = NaN(2,np);
            for p = 1:np
                
                % Find the last (perfect) integration parameter
                if np == 1, p = numel(A.Spec.ParamGrid{1}); end
                
                % Get the corresponding (scaled) tree
                pred = A.Fit.ScSurp{mmidx}(:,:,:,p,q,w);
                
                % Back to the surprise scale (it was saved in the scaled format)
                beta1 = A.Fit.Beta{mmidx}(1,p,q,w);
                beta0 = A.Fit.Beta{mmidx}(2,p,q,w);
                pred = (pred - beta0) ./ beta1;
                
                % Keep only the longest patterns
                pred = squeeze(pred(end,:,:));
                pred = reshape(pred, [n,1]);
                
                % Remove the "leave-one-out" observation
                pred = pred(idx);
                
                % Try to fit that
                warning('off');
                [MSE(p), ~, R2(p), ~, ~, ~, beta(:,p)] = Squires1976_FitModel(pred, loo);
                warning('on');
            end
            
            % Get the best parameter
            if     np == 1, bpi = numel(A.Spec.ParamGrid{1});
            elseif np > 1, [~,bpi] = min(MSE); end
            bestp(m,mm,j) = A.Spec.ParamGrid{mmidx}(bpi);
            
            % Get the corresponding (scaled) tree
            pred = A.Fit.ScSurp{mmidx}(:,:,:,bpi,q,w);
            
            % Keep only the longest patterns
            pred = squeeze(pred(end,:,:));
            pred = reshape(pred, [n,1]);

            % Get the prediction for the observation that was left out and
            % measure its distance to the data
            epsilon(m,mm,j) = Data(j) - pred(j);
            
            % Measure the variance explained between the current model
            % predictions (with its based parameter based on the loo fit)
            % and the predictions of the models on which the loo loop is
            % looping on
            R2within(m,mm,j) = R2(bpi);
            
            % Measure the variance explained between the current model
            % predictions (with its based parameter based on the loo fit)
            % and the Squires' data
            [~, ~, R2generalized(m,mm,j), ~, ~, ~, ~, ~] = Squires1976_FitModel(pred, Data);
        end
    end
    fprintf('\n');
end

%% PLOT THE GOODNESS OF FIT
%  ========================

% Choose which R2 to display
R2 = R2generalized;
%R2 = R2within;

% Indices of models to look at
perfidx = 1:3;
leakidx = 7:9;

% Average of R2 in specific models
R2avg_AtoB{1} = mean(R2(leakidx,leakidx),3); % leak to leak
R2avg_AtoB{2} = mean(R2(leakidx,perfidx),3); % leak to perfect
R2avg_AtoB{3} = mean(R2(perfidx,leakidx),3); % perfect to leak
R2avg_AtoB{4} = mean(R2(perfidx,perfidx),3); % perfect to perfect

% Standard deviation of R2 in specific models
R2std_AtoB{1} = std(R2(leakidx,leakidx),[],3); % leak to leak
R2std_AtoB{2} = std(R2(leakidx,perfidx),[],3); % leak to perfect
R2std_AtoB{3} = std(R2(perfidx,leakidx),[],3); % perfect to leak
R2std_AtoB{4} = std(R2(perfidx,perfidx),[],3); % perfect to perfect

% Average of best parameters in specific models
BPavg_AtoB{1} = mean(bestp(leakidx,leakidx),3); % leak to leak
BPavg_AtoB{2} = mean(bestp(leakidx,perfidx),3); % leak to perfect
BPavg_AtoB{3} = mean(bestp(perfidx,leakidx),3); % perfect to leak
BPavg_AtoB{4} = mean(bestp(perfidx,perfidx),3); % perfect to perfect

% Standard deviation of best parameters in specific models
BPstd_AtoB{1} = std(bestp(leakidx,leakidx),[],3); % leak to leak
BPstd_AtoB{2} = std(bestp(leakidx,perfidx),[],3); % leak to perfect
BPstd_AtoB{3} = std(bestp(perfidx,leakidx),[],3); % perfect to leak
BPstd_AtoB{4} = std(bestp(perfidx,perfidx),[],3); % perfect to perfect

from = {'leaky', 'leaky', 'perfect', 'perfect'};
to   = {'leaky', 'perfect', 'leaky', 'perfect'};

figure('Position', [0.1859 0.1683 0.6286 0.6250]);

pos = [-2/9, 0, 2/9];

for i = 1:4
    subplot(2,2,i);
    
    % Plot the bars
    bar(1:numel(stat), R2avg_AtoB{i}, 1, 'LineWidth', 1); hold('on');
    for m = 1:3
        for mm = 1:3
            
            % Display best parameter
            x = m + pos(mm);
            text(x, 0.1, sprintf('%1.0f', BPavg_AtoB{i}(m,mm)));
            text(x, 0.05, sprintf('(%.2f)', BPstd_AtoB{i}(m,mm)), 'FontAngle', 'Italic');
            
            % Plot error bars
            M = R2avg_AtoB{i}(m,mm);
            S = R2std_AtoB{i}(m,mm);
            plot(repmat(x,1,2), [M-S, M+S], 'k-', 'LineWidth', 1);
        end
    end
    
    ylim([0,1.2]); set(gca, 'YTick', 0:0.2:1); colormap(lines(3));
    
    set(gca, 'XTick', 1:numel(stat), 'XTickLabel', stat, 'XTickLabelRotation', 0);
    set(gca, 'Box', 'Off', 'LineWidth', 1, 'FontSize', 14);
    lgd = legend(stat, 'Box', 'Off', 'Location', 'North', 'Orientation', 'Horizontal'); 
    
    xlabel(sprintf('From %s ...', from{i}));
    title(sprintf('... to %s ...', to{i}), 'FontWeight', 'Normal');
    ylabel('< max_{params} R^{2} >_{loo} (+/- std)');
end
