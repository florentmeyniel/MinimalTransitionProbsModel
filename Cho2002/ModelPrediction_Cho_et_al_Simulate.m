% This script simulates a long, purely random sequence. It computes the
% posterior inference for different observer models and returns surprise 
% levels for each stimulus in the sequence and each model. We consider
% three models:
%   a learning of transition probabilities
%   a learning of stimulus frequency
%   a learning of alternation frequency
% All models have "fixed belief" with leaky integration.
% Surprise levels are then sorted based on patterns of five stimuli. These
% patterns can form, in total, 2^4=16 sequences of alternations and
% repetitions.
% The results are saved and should be compared with Cho et al. (2002) data.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% =========================================================================
%                               INITIALIZATION
% =========================================================================

clear; close all;
addpath('../IdealObserversCode/')
addpath('../Tools/')

L       = 1e3;              % length of the sequence of observation
N       = 4;                % length of the patterns

% Get label order from Cho et al
patterns_Cho = {
    'RRRR'
    'ARRR'
    'RARR'
    'AARR'
    'RRAR'
    'ARAR'
    'RAAR'
    'AAAR'
    'RRRA'
    'ARRA'
    'RARA'
    'AARA'
    'RRAA'
    'ARAA'
    'RAAA'
    'AAAA'};

% =========================================================================
%                               RUN SIMULATION
% =========================================================================

% simulate a purely random sequence
s = (rand(1,L+1)>0.5) + 1;

% compute the Ideal Observer inference based on transition probabilities
in.s                = s;                % sequence
in.learned          = 'transition';     % estimate transition
in.jump             = 0;                % estimate without jumps
in.opt.MemParam     = {'Decay', 16};    % memory limit
in.opt.AboutFirst   = 'WithoutFirst';   % discard 1st observation for analytical solutions
in.opt.ReturnDist   = 0;                % do not return full posterior distributions
in.opt.priorp1g2    = [1 1];            % uniform Beta prior
in.opt.priorp2g1    = [1 1];            % uniform Beta prior
in.verbose          = 1;                % check that no default values are used.
out{1} = IdealObserver(in);

in.s                = s;                % sequence
in.learned          = 'frequency';      % estimate transition
in.jump             = 0;                % estimate without jumps
in.opt.MemParam     = {'Decay', 10};    % memory limit
in.opt.ReturnDist   = 0;                % do not return full posterior distributions
in.opt.priorp1      = [1 1];            % uniform Beta prior
in.verbose          = 1;                % check that no default values are used.
out{2} = IdealObserver(in);

in.s                = abs(diff([1 s]))+1; % sequence translated into repeat (1) or alternate (2)
in.learned          = 'frequency';      % estimate transition
in.jump             = 0;                % estimate without jumps
in.opt.MemParam     = {'Decay', 25};    % memory limit
in.opt.ReturnDist   = 0;                % do not return full posterior distributions
in.opt.priorp1      = [1 1];            % uniform Beta prior
in.verbose          = 1;                % check that no default values are used.
out{3} = IdealObserver(in);


% =========================================================================
%                           SORT PATTERNS
% =========================================================================

% convert the sequence into repeat (1/A) and alternate (2/B)
% note that 1/2 correspond to labels in the sequence, and A/B to the labels
% as returned by GetPatLL
sRep = abs(diff(s))+1;

for iObs = 1:3
    % search likelihood of repeat in patterns of N
    [mLL_R, ~, ~, ~, ~, ~, ~, seqnames_R] = ...
        GetPatLL(N, sRep, out{iObs}.surprise(2:end), false, 1, 0);
    
    % search likelihood of alternate  in patterns of N
    [mLL_A, ~, ~, ~, ~, ~, ~, seqnames_A] = ...
        GetPatLL(N, sRep, out{iObs}.surprise(2:end), false, 2, 0);
    
    % concatenate results in a single cell
    rescell = cell((2^(N-1))*2, N+1);
    rescell(1:2^(N-1), 1:N) = seqnames_R{N};
    rescell(2^(N-1)+1:(2^(N-1))*2, 1:N) = seqnames_A{N};
    rescell(1:2^(N-1), N+1) = num2cell(mLL_R(end,:)', 2);
    rescell(2^(N-1)+1:(2^(N-1))*2, N+1) = num2cell(mLL_A(end,:)', 2);
    
    % get pattern labels from this cell
    pattern_rescell = cell((2^(N-1))*2, 1);
    for iPat = 1:16
        pattern_rescell{iPat} = [rescell{iPat,1}, rescell{iPat,2}, rescell{iPat,3}, rescell{iPat,4}];
    end
    
    % =========================================================================
    %                           PLOT RESULTS
    % =========================================================================
    
    % Convert the pattern in pattern_rescell (A/B) into Cho et al format (R/A)
    patterns_GetPatLL = cell(16,1);
    for iPat = 1:16
        patterns_GetPatLL{iPat,1} = 'xxxx';
        for iLetter = 1:4
            if strcmp(patterns_Cho{iPat,1}(iLetter), 'R'); patterns_GetPatLL{iPat,1}(iLetter) = 'A'; end
            if strcmp(patterns_Cho{iPat,1}(iLetter), 'A'); patterns_GetPatLL{iPat,1}(iLetter) = 'B'; end
        end
    end
    
    % get likelihood for each pattern (same order as in Cho et al)
    IOprediction_Cho_order = zeros(1, 16);
    for iPat = 1:16
        ind = strncmp(patterns_GetPatLL{iPat}, pattern_rescell, 5);
        IOprediction_Cho_order(iPat) = rescell{ind, end};
    end
    
    % store data
    if iObs == 1; Prediction_TransitionProbs =  IOprediction_Cho_order; end
    if iObs == 2; Prediction_StimProb =  IOprediction_Cho_order; end
    if iObs == 3; Prediction_AltProb =  IOprediction_Cho_order; end    
end

% save data
timestamp = sprintf('%d-%d-%d_%d-%d-%1.0f_%5.0f', clock, rand*1e5);
save(sprintf('Simulation_Cho_%s.mat', timestamp), ...
    'Prediction_TransitionProbs', 'Prediction_StimProb', 'Prediction_AltProb', ...
    'L', 'N', 'patterns_Cho')
