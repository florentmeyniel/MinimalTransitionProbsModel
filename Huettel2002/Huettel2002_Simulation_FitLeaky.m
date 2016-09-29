% This script simulates a long random sequence, corresponding to purely
% unpredictable binary outcomes. It is like a series of fair coin flips.
% The script estimates the prediction of three ideal observers that learn,
% respectively, the frequency of outcomes, the frequency of alternation
% between outcomes and the transition probabilities between outcomes.
% All observers have leaky integration and learn mostly from recent
% observations. Different leak values are tested for data fitting.
% Then, we look for local series of repeated and alternated elements, and
% there violation. A series of repeated elements is violated if a
% new observation corresponds to a change, and a series of alternating
% elements is violated if a new observation corresponds to a repetition.
% The ideal observers' surprise levels are sorted in term of continuation
% and violation of these local series of repetitions and alternations.
%
% This theoretical analysis should be compared to the behavioral and
% fMRI results by Huettel et al, Nature Neuroscience 2002.
% The result of this simulation are saved.
%
% Note that this code loops over 3 observer types, and for each observer,
% over many values of the leak parameter. Two modes are proposed: parallel
% computation over observer type, or serial computation. The parallel mode
% may not work on all operating systems, and it requires to specify the 
% command line to invoke matlab from a terminal.
%
% Copyright 2016 Florent Meyniel & Maxime Maheu

% Clear the place
clear; close('all');
addpath ../IdealObserversCode/
addpath ../Tools/

% compute in parallel (=1) or not (=0)?
do_parallel = 1;

% Generate a long, maximally unpredictable binary sequence (= series of
% coin flips)
p1  = 0.5;
L   = 1e5;
s   = GenRandSeq(L, p1);

% parameter grid
grid_leak = 1:25;

% Common options
in.jump             = 0;                % assume that there are no "jumps" (changes in generative statistics)
in.opt.AboutFirst   = 'WithoutFirst';   % discard 1st observation for analytical solutions
n                   = 20;               % resolution of the univariate probability grid
in.opt.pgrid        = linspace(0,1,n);  % grid to return full distributions
in.opt.ReturnDist   = 0;                % Do not return full posterior distributions
in.opt.priorp1      = [1 1];            % uniform Beta prior on frequency
in.opt.priorp1g2    = [1 1];            % uniform Beta prior on p(1|2)
in.opt.priorp2g1    = [1 1];            % uniform Beta prior on p(2|1)
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
            fname       = 'LeakItemFreq_param.mat';
        elseif iObs == 3
            in.learned  = 'frequency';
            in.s        = [1 abs(diff(s))+1];   % recode sequence as alternation / repetition
            fname       = 'LeakyAltFreq_param.mat';
        end
        save(fname, 'in', 'grid_leak')
        
        cmd = [sprintf('%s -nosplash -nodesktop -r "', 'matlab'), ...        % open Matlab
            sprintf('tic; '), ...                                                   % tic
            sprintf('cd %s ; ', pwd), ...                                           % go into current directory
            sprintf(['try Huettel2002_Simulation_1Obs_MultParam(''%s'');' ...
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
    JobList = {'LeakyTransition_param', 'LeakItemFreq_param', 'LeakyAltFreq_param'};
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
            fname       = 'LeakItemFreq_param.mat';
        elseif iObs == 3
            in.learned  = 'frequency';
            in.s        = [1 abs(diff(s))+1];   % recode sequence as alternation / repetition
            fname       = 'LeakyAltFreq_param.mat';
        end
        save(fname, 'in', 'grid_leak')
        
        Huettel2002_Simulation_1Obs_MultParam(fname)
    end
end

ObsName             = ...
    {'transition probs.', ...
    'stimulus freq.', ...
    'alternation freq.', ...
    };

% =========================================================================
%                 Huettel et al. Nat Neuro 2002
% =========================================================================

% nomenclature of Huettel:
% sequence length = 0 => mean for all events
% sequence length = 1 => x1 [x2] y      (with x1 # x2)
% sequence length = 2 => x1 [x2 x3] y   (with x1 # x2)
% the data correspond what happens at event y

% Compute transitions between events
trn = diff(s);
maxrep = 8;
n_param = numel(grid_leak);

% get mean estimates for each observer
mrep_cont = nan(3, maxrep, n_param);
mrep_viol = nan(3, maxrep, n_param);
srep_cont = nan(3, maxrep, n_param);
srep_viol = nan(3, maxrep, n_param);
malt_cont = nan(3, maxrep, n_param);
malt_viol = nan(3, maxrep, n_param);
salt_cont = nan(3, maxrep, n_param);
salt_viol = nan(3, maxrep, n_param);

% load the simulated data
tmp = load('LeakyTransition_res.mat');  ObsLL{1} = tmp.ObsLL;
tmp = load('LeakItemFreq_res.mat');     ObsLL{2} = tmp.ObsLL;
tmp = load('LeakyAltFreq_res.mat');     ObsLL{3} = tmp.ObsLL;

% 1: LOOK FOR REPETITIONS
% =======================
maxrep = 8;
repA_viol = cell(3,maxrep); repB_viol = cell(3,maxrep); rep_viol = cell(3,maxrep);
repA_cont = cell(3,maxrep); repB_cont = cell(3,maxrep); rep_cont = cell(3,maxrep);
for k = 1:maxrep
    for iObs = 1:3
        repA_viol{iObs,k} = [];
        repB_viol{iObs,k} = [];
        rep_viol{iObs,k}  = [];
        repA_cont{iObs,k} = [];
        repB_cont{iObs,k} = [];
        rep_cont{iObs,k}  = [];
    end
end

% VIOLATION OF REPETITION
% ~~~~~~~~~~~~~~~~~~~~~~~
% COMPUTE FOR LENGTH = 1 (i.e. Y [X] Y)
for k = 1:length(trn)-1
    % look for [1 -1] transition (A_{k} B_{k+1} A_{k+2})
    if trn(k) == 1 && trn(k+1) == -1
        for iObs = 1:3
            repA_viol{iObs,1} = [repA_viol{iObs,1}, ObsLL{iObs}(:,k+2)];
        end
    end
    
    % look for [-1 1] transition (B_{k} A_{k+1} B_{k+2})
    if trn(k) == -1 && trn(k+1) == 1
        for iObs = 1:3
            repB_viol{iObs,1} = [repB_viol{iObs,1}, ObsLL{iObs}(:,k+2)];
        end
    end
end
for iObs = 1:3
    rep_viol{iObs,1} = [repA_viol{iObs,1}, repB_viol{iObs,1}];
end


% COMPUTE FOR LENGTH > 1 (i.e. Y [X ... X] Y)
% get the rank of the 1st element of the first repetitions
k = find(trn == 0, 1, 'first');
while k < length(trn) - 1
    lrep = find(trn(k:end), 1, 'first'); % lrep-1 repetition
    k = k+lrep;                          % index of violation after repetitions
    
    if lrep <= maxrep
        for iObs = 1:3
            if s(k) == 1, repA_viol{iObs,lrep} = [repA_viol{iObs,lrep}, ObsLL{iObs}(:,k)];     % A violated a local repetition
            else          repB_viol{iObs,lrep} = [repB_viol{iObs,lrep}, ObsLL{iObs}(:,k)]; end % B violated a local repetition
            rep_viol{iObs,lrep} = [rep_viol{iObs,lrep}, ObsLL{iObs}(:,k)];
        end
    end
    
    % from then, get the index of the 1st element of a repetition
    k = k + find(trn(k:end) == 0, 1, 'first') - 1;
end

% CONTINUATION OF REPETITION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
% get the rank of the 1st element of the first repetition
k = find(trn == 0, 1, 'first');
while k < length(trn) - 1
    lrep = find(trn(k:end), 1, 'first'); % if lrep=1: there is a change; if lrep>1: there are lrep-1 repetitions
    
    % Get surprise for each repetition within this repetition chunk
    for kk = 1:min(lrep-1, maxrep)
        for iObs = 1:3
            if s(k+kk) == 1, repA_cont{iObs,kk} = [repA_cont{iObs,kk}, ObsLL{iObs}(:,k+kk)];     % A repetition locally continued
            else             repB_cont{iObs,kk} = [repB_cont{iObs,kk}, ObsLL{iObs}(:,k+kk)]; end % B repetition locally continued
            rep_cont{iObs,kk} = [rep_cont{iObs,kk}, ObsLL{iObs}(:,k+kk)];
        end
    end
    
    % go to the 1st element after this chunk (i.e. index of first
    % change)
    k = k + lrep;
    
    % from then, get the index of the 1st element of a repetition
    k = k + find(trn(k:end) == 0, 1, 'first') - 1;
end

% Compute mean surprise levels
for k = 1:maxrep
    for iObs = 1:3
        if ~isempty(rep_cont{iObs,k})
            mrep_cont(iObs, k, :) = mean(rep_cont{iObs,k},2);
            srep_cont(iObs, k, :) = stderror(rep_cont{iObs,k},2);
        end
        if ~isempty(rep_viol{iObs,k})
            mrep_viol(iObs, k, :) = mean(rep_viol{iObs,k},2);
            srep_viol(iObs, k, :) = stderror(rep_viol{iObs,k},2);
        end
    end
end

% 2: LOOK FOR ALTERNATION
% ====================
altA_viol = cell(3,maxrep); altB_viol = cell(3,maxrep); alt_viol = cell(3,maxrep);
altA_cont = cell(3,maxrep); altB_cont = cell(3,maxrep); alt_cont = cell(3,maxrep);
for k = 1:maxrep
    for iObs = 1:3
        altA_viol{iObs,k} = [];
        altB_viol{iObs,k} = [];
        alt_viol{iObs,k}  = [];
        altA_cont{iObs,k} = [];
        altB_cont{iObs,k} = [];
        alt_cont{iObs,k}  = [];
    end
end

% CONTINUATION OF ALTERNATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% get the rank of the 1st non repeated element
k = find(trn ~= 0, 1, 'first');
while k < length(trn) - 1
    lnrep = find(trn(k:end) == 0, 1, 'first'); % lnrep-1 alternations
    
    for kk = 1:min(lnrep-1, maxrep)
        for iObs = 1:3
            if s(k+kk) == 1, altA_cont{iObs,kk} = [altA_cont{iObs,kk}, ObsLL{iObs}(:,k+kk)];     % A violated a local alternation
            else             altB_cont{iObs,kk} = [altB_cont{iObs,kk}, ObsLL{iObs}(:,k+kk)]; end % B violated a local alternation
            alt_cont{iObs,kk} = [alt_cont{iObs,kk}, ObsLL{iObs}(:,k+kk)];
        end
    end
    
    % go to the 1st element after this chunk (i.e. to the next first
    % repeated event)
    k = k + lnrep;
    
    % from then, get the rank of the 1st element of a alternation
    k = k + find(trn(k:end) ~= 0, 1, 'first') - 1;
end

% VIOLATION OF ALTERNATION
% ~~~~~~~~~~~~~~~~~~~~~~~~
% COMPUTE FOR LENGTH = 1 (i.e. X [X] X)
for k = 1:length(s) - 2
    if trn(k) == 0 && trn(k+1) == 0
        for iObs = 1:3
            if s(k+2) == 1, altA_viol{iObs,1} = [altA_viol{iObs,1}, ObsLL{iObs}(:,k+2)];
            else            altB_viol{iObs,1} = [altB_viol{iObs,1}, ObsLL{iObs}(:,k+2)]; end
            alt_viol{1} = [alt_viol{iObs,1}, ObsLL{iObs}(:,k+2)];
        end
    end
end

% COMPUTE FOR LENGTH > 1
% get the rank of the 1st non repeated element
k = find(trn ~= 0, 1, 'first');
while k < length(trn) - 1
    lnrep = find(trn(k:end) == 0, 1, 'first'); % lnrep-1 alternetions
    
    % go to the 1st element after this chunk of repetitions
    k = k + lnrep;
    
    if lnrep <= maxrep
        for iObs = 1:3
            if s(k) == 1, altA_viol{iObs,lnrep} = [altA_viol{lnrep}, ObsLL{iObs}(:,k)];     % A violated a local alternation
            else          altB_viol{iObs,lnrep} = [altB_viol{lnrep}, ObsLL{iObs}(:,k)]; end % B violated a local alternation
            alt_viol{iObs,lnrep} = [alt_viol{iObs,lnrep}, ObsLL{iObs}(:,k)];
        end
    end
    
    % from then, get the index of the 1st element of an alternation
    k = k + find(trn(k:end) ~= 0, 1, 'first') - 1;
end

% Compute mean surprise levels
for k = 1:maxrep
    for iObs = 1:3
        if ~isempty(alt_cont{iObs,k})
            malt_cont(iObs, k, :) = mean(alt_cont{iObs,k},2);
            salt_cont(iObs, k, :) = stderror(alt_cont{iObs,k},2);
        end
        if ~isempty(alt_viol{iObs,k})
            malt_viol(iObs, k, :) = mean(alt_viol{iObs,k},2);
            salt_viol(iObs, k, :) = stderror(alt_viol{iObs,k},2);
        end
    end
end


timestamp = sprintf('%d-%d-%d_%d-%d-%1.0f_%5.0f', clock, rand*1e5);
save(sprintf('Huettel2002_Simulation_Results_%s.mat', timestamp), ...
    'mrep_cont', 'srep_cont', ...
    'mrep_viol', 'srep_viol', ...
    'malt_cont', 'salt_cont', ...
    'malt_viol', 'salt_viol', ...
    'maxrep', 'ObsName', 'grid_leak')
