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

L       = 1e5;              % length of the sequence of observation
N       = 4;                % length of the patterns

% parameter grid
grid_leak = 1:25;

% compute in parallel (=1) or not (=0)?
do_parallel = 1;

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

% Common options
in.jump             = 0;                % estimate without jumps
in.opt.AboutFirst   = 'WithoutFirst';   % discard 1st observation for analytical solutions
in.opt.ReturnDist   = 0;                % do not return full posterior distributions
in.opt.priorp1      = [1 1];            % uniform Beta prior on frequency
in.opt.priorp1g2    = [1 1];            % uniform Beta prior
in.opt.priorp2g1    = [1 1];            % uniform Beta prior
in.verbose          = 1;                % check that no default values are used.

% RUN THE IDEAL OBSERVER ESTIMATION
switch do_parallel
    case 1
        % Run the estimation in parallel
        timestamp = sprintf('%d-%d-%d_%d-%d-%1.0f_%5.0f', clock, rand*1e5);
        for iObs = 1:3
            if iObs == 1
                in.learned  = 'transition';
                in.s        = s;
                fname       = 'LeakyTransition_param.mat';
            elseif iObs == 2
                in.learned  = 'frequency';
                in.s        = s;
                fname       = 'LeakyItemFreq_param.mat';
            elseif iObs == 3
                in.learned  = 'frequency';
                in.s        = [1 abs(diff(s))+1];   % recode sequence as alternation / repetition
                fname       = 'LeakyAltFreq_param.mat';
            end
            save(fname, 'in', 'grid_leak')
            
            cmd = [sprintf('%s -nosplash -nodesktop -r "', 'matlab'), ...        % open Matlab
                sprintf('tic; '), ...                                                   % tic
                sprintf('cd %s ; ', pwd), ...                                           % go into current directory
                sprintf(' addpath ../Tools/ ; '), ...                                   % set path
                sprintf(['try Simulation_1Obs_MultParam(''%s'');' ...
                ' catch ME; disp(getReport(ME,''extended'')); end ; '], ...             % execute function with job
                fname) ...
                sprintf('mkdir(''over_%s_%s'') ; ', timestamp, fname(1:end-4)), ...  % make a temporary directory to signal job end
                sprintf('toc; '), ...                                                   % toc
                'quit' ...                                                              % close Matlab
                '" '...
                sprintf('> log_file_%s.txt ', fname(1:end-4)), ...
                '&'];
            unix(cmd)
        end
        
        % Update the list of batches currently running
        JobList = {'LeakyTransition_param', 'LeakyItemFreq_param', 'LeakyAltFreq_param'};
        RunningJobs = 1:3;
        while ~isempty(RunningJobs)
            for iJob = RunningJobs
                % check if job has returned (with the 'probe directory')
                if exist( sprintf('over_%s_%s', timestamp, JobList{iJob}), 'dir')
                    % remove the probe
                    rmdir(sprintf('over_%s_%s', timestamp, JobList{iJob}))
                    
                    % Remove job from list of running jobs
                    RunningJobs = setdiff(RunningJobs, iJob);
                    fprintf('\n Job completed: %s', JobList{iJob})
                end
            end
        end
        
    case 0
        % Run the observer one at a time
        for iObs = 1:3
            if iObs == 1
                in.learned  = 'transition';
                in.s        = s;
                fname       = 'LeakyTransition_param.mat';
            elseif iObs == 2
                in.learned  = 'frequency';
                in.s        = s;
                fname       = 'LeakyItemFreq_param.mat';
            elseif iObs == 3
                in.learned  = 'frequency';
                in.s        = [1 abs(diff(s))+1];   % recode sequence as alternation / repetition
                fname       = 'LeakyAltFreq_param.mat';
            end
            save(fname, 'in', 'grid_leak')
            
            Simulation_1Obs_MultParam(fname)
        end
end

% =========================================================================
%                           SORT PATTERNS
% =========================================================================

% convert the sequence into repeat (1/A) and alternate (2/B)
% note that 1/2 correspond to labels in the sequence, and A/B to the labels
% as returned by GetPatLL
sRep = abs(diff(s))+1;

% load the simulated data
tmp = load('LeakyTransition_res.mat');  ObsLL{1} = tmp.ObsLL;
tmp = load('LeakyItemFreq_res.mat');     ObsLL{2} = tmp.ObsLL;
tmp = load('LeakyAltFreq_res.mat');     ObsLL{3} = tmp.ObsLL;

% Search repeat patterns
[ ~, ~, ~, ~, ~, getLLpos_R, ~, seqnames_R] = ...
    GetPatLL( N, sRep, zeros(1,numel(sRep)), false, 1, 0);

% Search alternate patterns
[ ~, ~, ~, ~, ~, getLLpos_A, ~, seqnames_A] = ...
    GetPatLL(N, sRep, zeros(1,numel(sRep)), false, 2, 0);

% number of pattern for maximal lenght
maxlen = size(ff2n(N),1)/2;
for iObs = 1:3
    
    IOprediction_Cho_order = nan(numel(grid_leak,16));
    fprintf('\n Searching patterns for observer %d/%d', iObs, 3)
    
    for iParam = 1:numel(grid_leak)
        fprintf('\n   param %d/%d', iParam, numel(grid_leak))
        
        % log-likelihood for this observer and this parameter
        seqLL = ObsLL{iObs}(iParam,2:end);
        
        mLL_R = NaN(N,maxlen);
        mLL_A = NaN(N,maxlen);
        for kk = 1:N
            for k = 1:2^(kk-1)
                % Mean likelihood of the pattern
                mLL_R(kk,k) = nanmean(seqLL(getLLpos_R{kk,k})); 
                mLL_A(kk,k) = nanmean(seqLL(getLLpos_A{kk,k})); 
            end
        end
        
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
        for iPat = 1:16
            ind = strncmp(patterns_GetPatLL{iPat}, pattern_rescell, 5);
            IOprediction_Cho_order(iParam,iPat) = rescell{ind, end};
        end
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
    'L', 'N', 'patterns_Cho', 'grid_leak')
