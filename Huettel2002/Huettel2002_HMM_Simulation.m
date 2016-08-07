% This script simulates a long random sequence, corresponding to purely
% unpredictable binary outcomes. It is like a series of fair coin flips.
% The script estimates the prediction of three ideal observers that learn, 
% respectively, the frequency of outcomes, the frequency of alternation 
% between outcomes and the transition probabilities between outcomes. 
% Both observers assume that the statistic they infer
% may change over time, so that they learn mostly from recent observations.
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
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% Clear the place
clear; close('all');
addpath ../IdealObserversCode/
addpath ../Tools/


% Generate a long, maximally unpredictable binary sequence (= series of
% coin flips)
p1  = 0.5;
L   = 1e5; 
s   = GenRandSeq(L, p1);

% Common options
in.s                = s;                % sequence
in.jump             = 1;                % assume that there are "jumps" (changes in generative statistics)
in.mode             = 'HMM';            % compute with HMM
n                   = 20;               % resolution of the univariate probability grid
in.opt.pgrid        = linspace(0,1,n);  % grid to return full distributions
in.verbose          = 1;                % check that no default values are used.

% Run the Markov observer (learn transition probabilities between observations)
in.learned          = 'transition';
in.opt.Alpha0       = ones(n)/(n^2);    % uniform prior on transition probabilities
in.opt.pJ           = 0.013;            % prior probability of jump
out                 = IdealObserver(in);
ObsLL{1}            = out.surprise;

% Run the Bernoulli observer (learn frequency of observation)
in.learned          = 'frequency';
in.opt.Alpha0       = ones(1,n)/n;      % uniform prior on frequency
in.opt.pJ           = 0.019;            % prior probability of jump
out                 = IdealObserver(in);
ObsLL{2}            = out.surprise;

% Run the Bernoulli observer that learns the frequency of transitions
in.s                = [1 abs(diff(s))+1];   % recode sequence as alternation / repetition
in.learned          = 'frequency';
in.opt.Alpha0       = ones(1,n)/n;      % uniform prior on frequency
in.opt.pJ           = 0.019;            % prior probability of jump
out                 = IdealObserver(in);
ObsLL{3}            = out.surprise;

ObsName             = {'transition probs.', ...
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

% get mean estimates for each observer
mrep_cont = nan(3,maxrep);
mrep_viol = nan(3,maxrep);
srep_cont = nan(3,maxrep);
srep_viol = nan(3,maxrep);
malt_cont = nan(3,maxrep);
malt_viol = nan(3,maxrep);
salt_cont = nan(3,maxrep);
salt_viol = nan(3,maxrep);
for iObs = 1:3
    
    % Get surprise for this observer
    LL = ObsLL{iObs};
    
    % 1: LOOK FOR REPETITIONS
    % =======================
    maxrep = 8;
    repA_viol = cell(maxrep,1); repB_viol = cell(maxrep,1); rep_viol = cell(maxrep,1);
    repA_cont = cell(maxrep,1); repB_cont = cell(maxrep,1); rep_cont = cell(maxrep,1);
    for k = 1:maxrep
        repA_viol{k} = [];
        repB_viol{k} = [];
        rep_viol{k}  = [];
        repA_cont{k} = [];
        repB_cont{k} = [];
        rep_cont{k}  = [];
    end
    
    % VIOLATION OF REPETITION
    % ~~~~~~~~~~~~~~~~~~~~~~~
    % COMPUTE FOR LENGTH = 1 (i.e. Y [X] Y)
    for k = 1:length(trn)-1
        % look for [1 -1] transition (A_{k} B_{k+1} A_{k+2})
        if trn(k) == 1 && trn(k+1) == -1
            repA_viol{1} = [repA_viol{1}, LL(k+2)];
        end
        
        % look for [-1 1] transition (B_{k} A_{k+1} B_{k+2})
        if trn(k) == -1 && trn(k+1) == 1
            repB_viol{1} = [repB_viol{1}, LL(k+2)];
        end
    end
    rep_viol{1} = [repA_viol{1}, repB_viol{1}];
    
    
    % COMPUTE FOR LENGTH > 1 (i.e. Y [X ... X] Y)
    % get the rank of the 1st element of the first repetitions
    k = find(trn == 0, 1, 'first');
    while k < length(trn) - 1
        lrep = find(trn(k:end), 1, 'first'); % lrep-1 repetition
        k = k+lrep;                          % index of violation after repetitions
        
        if lrep <= maxrep
            if s(k) == 1, repA_viol{lrep} = [repA_viol{lrep}, LL(k)];     % A violated a local repetition
            else          repB_viol{lrep} = [repB_viol{lrep}, LL(k)]; end % B violated a local repetition
            rep_viol{lrep} = [rep_viol{lrep}, LL(k)];
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
            if s(k+kk) == 1, repA_cont{kk} = [repA_cont{kk}, LL(k+kk)];     % A repetition locally continued
            else             repB_cont{kk} = [repB_cont{kk}, LL(k+kk)]; end % B repetition locally continued
            rep_cont{kk} = [rep_cont{kk}, LL(k+kk)];
        end
        
        % go to the 1st element after this chunk (i.e. index of first
        % change)
        k = k + lrep;
        
        % from then, get the index of the 1st element of a repetition
        k = k + find(trn(k:end) == 0, 1, 'first') - 1;
    end
    
    % Compute mean surprise levels
    for k = 1:maxrep
        if ~isempty(rep_cont{k})
            mrep_cont(iObs, k) = mean(rep_cont{k});
            srep_cont(iObs, k) = stderror(rep_cont{k});
        end
        if ~isempty(rep_viol{k})
            mrep_viol(iObs, k) = mean(rep_viol{k});
            srep_viol(iObs, k) = stderror(rep_viol{k});
        end
    end
    
    % 2: LOOK FOR ALTERNATION
    % ====================
    altA_viol = cell(maxrep,1); altB_viol = cell(maxrep,1); alt_viol = cell(maxrep,1);
    altA_cont = cell(maxrep,1); altB_cont = cell(maxrep,1); alt_cont = cell(maxrep,1);
    for k = 1:maxrep
        altA_viol{k} = [];
        altB_viol{k} = [];
        alt_viol{k}  = [];
        altA_cont{k} = [];
        altB_cont{k} = [];
        alt_cont{k}  = [];
    end
    
    % CONTINUATION OF ALTERNATION
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % get the rank of the 1st non repeated element
    k = find(trn ~= 0, 1, 'first');
    while k < length(trn) - 1
        lnrep = find(trn(k:end) == 0, 1, 'first'); % lnrep-1 alternations 
        
        for kk = 1:min(lnrep-1, maxrep)
            if s(k+kk) == 1, altA_cont{kk} = [altA_cont{kk}, LL(k+kk)];     % A violated a local alternation
            else             altB_cont{kk} = [altB_cont{kk}, LL(k+kk)]; end % B violated a local alternation
            alt_cont{kk} = [alt_cont{kk}, LL(k+kk)];
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
            if s(k+2) == 1, altA_viol{1} = [altA_viol{1}, LL(k+2)];
            else            altB_viol{1} = [altB_viol{1}, LL(k+2)]; end
            alt_viol{1} = [alt_viol{1}, LL(k+2)];
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
            if s(k) == 1, altA_viol{lnrep} = [altA_viol{lnrep}, LL(k)];     % A violated a local alternation
            else          altB_viol{lnrep} = [altB_viol{lnrep}, LL(k)]; end % B violated a local alternation
            alt_viol{lnrep} = [alt_viol{lnrep}, LL(k)];
        end
        
        % from then, get the index of the 1st element of an alternation
        k = k + find(trn(k:end) ~= 0, 1, 'first') - 1;
    end
    
    % Compute mean surprise levels
    for k = 1:maxrep
        if ~isempty(alt_cont{k})
            malt_cont(iObs, k) = mean(alt_cont{k});
            salt_cont(iObs, k) = stderror(alt_cont{k});
        end
        if ~isempty(alt_viol{k})
            malt_viol(iObs, k) = mean(alt_viol{k});
            salt_viol(iObs, k) = stderror(alt_viol{k});
        end
    end
    
end

timestamp = sprintf('%d-%d-%d_%d-%d-%1.0f_%5.0f', clock, rand*1e5);
save(sprintf('Huettel2002_SimulationHMM_Results_%s.mat', timestamp), ...
    'mrep_cont', 'srep_cont', ...
    'mrep_viol', 'srep_viol', ...
    'malt_cont', 'salt_cont', ...
    'malt_viol', 'salt_viol', ...
    'maxrep', 'ObsName')
