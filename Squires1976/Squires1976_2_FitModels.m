% Fit the models to Squires et al 1976 data (or Kolossa et al 2013 data).
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
%% INITIALIZATION
%  ==============

% Clear the place
clear; close('all'); clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE SOME OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% On which dataset
d = 'SquiresScience1976';
%d = 'KolossaFIHN2013';

% Choose on which sequences to average (useful to check that simulations
% reached the asymptotic values)
seqidx = @(Nseqmax) 1:Nseqmax; % all sequences
%seqidx = @(Nseqmax) 1:2:Nseqmax; % only odd sequences
%seqidx = @(Nseqmax) 2:2:Nseqmax; % only even sequences

% Define the function to apply on weigths in the regression
%weightfun = @(x) sqrt(x); % \sqrt{p(pattern|p(1))}
weightfun = @(x) (x); % p(pattern|p(1))

% Define the patterns used to measure the quality of fit
pattofit = @(Npatmax) Npatmax; % only the longest patterns
%pattofit = @(Npatmax) 1:Npatmax; % all available patterns (i.e. the entire tree)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the functions
addpath(genpath('Squires1976'));
try cd('Squires1976'); catch, end

%% LOAD DATA
%  =========

% Load the results
load('Squires1976_SimulateModels.mat');

% Add informations
Info.DataSet   = d;
Info.SeqIdx    = seqidx;
Info.WeightFun = weightfun;
Info.FitPatLen = pattofit;

% Load data
[Data, nPat, p1] = Squires1976_PrepareData(sprintf('Data%s.mat', d));
nP = numel(p1);

% Reorder the blocks
[~,I] = sort(p1, 2, 'descend');
Data = Data(:,:,I);
p1 = p1(I);

% Get patterns' expected frequency
PatProb = NaN(size(Data));
for p = 1:nP
    PatProb(:,:,p) = GetPatProb(nPat, p1(p));
end

%% SCALE THE SIMULATIONS TO THE DATA
%  =================================

% Display what we are doing
fprintf('FINDING PATTERNS IN THE SEQUENCES AND SCALING THE AVERAGED SURPRISE TO THE DATA...\n\n');

% Get the number of estimated models
nMod = numel(Simu);

% Prepare variables for the quality of fit
Fit.MSE    = cell(1,nMod);
Fit.R2     = cell(1,nMod);
Fit.Beta   = cell(1,nMod);
Fit.ScSurp = cell(1,nMod);

% Custom regression options
nSeq = size(Spec.StimSeq,1);
seqidx = seqidx(nSeq);
pattofit = pattofit(nPat);

% For each model
for m = 1:nMod
    
    % Get the model
    stat = Info.Models{m,1};
    intg = Info.Models{m,2};
    fprintf('Learning "%s" with a "%s" integration:\n', stat, intg);

    % Prepare outputs
    [~, ~, nParam, nQuty, nSeq] = size(Simu{m});
    Fit.MSE{m}    = NaN(nParam, nQuty, 2);
    Fit.R2{m}     = NaN(nParam, nQuty, 2);
    Fit.Beta{m}   = NaN(2, nParam, nQuty, 2);
    Fit.ScSurp{m} = NaN([size(PatProb), nQuty, 2]);
    
    % Variables used for printing purpose
    np = numel(sprintf('%i', nParam));
    nd = numel(sprintf('%i', nSeq));
    rp = 11+(2*nd);
    bk = repmat('\b', [1,rp]);

    % For each quantity to estimate and each parameter
    for p = 1:nParam
        
        % Get the parameter value
        if ~isempty(Spec.ParamGrid{m}), param = Spec.ParamGrid{m}(p);
        else param = NaN; end
        fprintf([' - Parameter %', num2str(np), '.0f/%', num2str(np), ...
            '.0f (%g): \t%s'], p, nParam, param, repmat(' ', [1,rp]));
        
        % For each stimuli sequence, each global probability of 1 and
        % each quantity to estimate
        for q = 1:nQuty
            surp = NaN([size(PatProb), nSeq]);
            nobs = NaN([size(PatProb), nSeq]);
            for seq = 1:nSeq
                fprintf([bk, 'sequence %', num2str(nd), '.0f/%' , num2str(nd), '.0f\n'], seq, nSeq);
                for blk = 1:nP

                    % Look for specific patterns in the sequence
                    [surp(:,:,blk,seq), ...
                     nobs(:,:,blk,seq)] = ...
                        GetPatLL(nPat, ... % the length of patterns to look for
                                 Spec.StimSeq(seq,:,blk), ... % the stimuli sequence
                                 Simu{m}(:,blk,p,q,seq), ... % the simulation vector
                                 false, 1); % options
                end
            end
            
            % Keep only part of the data (useful for checking that the
            % simulations reached the asymptotic regime)
            surp = surp(:,:,:,seqidx);
            nobs = nobs(:,:,:,seqidx);
            
            % Average pattern-evoked surprise over sequences by doing
            % a weighted average whose weights are the number of
            % occurrence of each pattern in each sequence
            Surprise = nansum(surp .* nobs, 4) ./ nansum(nobs, 4);

            % Perform a linear regression between the data and the
            % surprise evoked by the longest available patterns only
            [Fit.MSE{m}(p,q,1), ...   % unweighted mean squared error
             Fit.MSE{m}(p,q,2), ...   % weighted mean squared error
             Fit.R2{m}(p,q,1), ...    % unweighted fraction of explained variance
             Fit.R2{m}(p,q,2), ...    % weighted fraction of explained variance
             ~, ~, ...
             Fit.Beta{m}(:,p,q,1), ...   % betas from the unweighted regression
             Fit.Beta{m}(:,p,q,2)] = ... % betas from the weighted regression
             Squires1976_FitModel(...
                Surprise(pattofit,:,:), ...        % the P300 amplitude evoked by each pattern
                Data(pattofit,:,:), ...            % the averaged surprise for each pattern
                weightfun(PatProb(pattofit,:,:))); % the probability of observing each pattern

            % Scale the other (shorter) patterns using the betas from
            % the previous regression
            for w = 1:2
                beta1 = Fit.Beta{m}(1,p,q,w);
                beta0 = Fit.Beta{m}(2,p,q,w);
                Fit.ScSurp{m}(:,:,:,p,q,w) = (Surprise .* beta1) + beta0;
            end
        end
    end
end

%% SAVE THE DATA
%  =============

save(sprintf('Squires1976_FitModels_%s.mat', d), '-v7.3', 'Info', 'Spec', 'Fit');
