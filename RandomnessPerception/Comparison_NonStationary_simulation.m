% This script simulates many sequences with different apparent 
% probabilities of alternation. The ideal observer inference is computed 
% assuming fixed belief and leaky integration. We distinguish between 
% learning transition probabilities, learning the stimulus frequency and 
% learning the frequency of alternations. We also compare different prior 
% belief on the frequency of changes (from rare to frequent). The entropy 
% of the last posterior prediction is computed to quantify the 
% unpredictability of the sequence of observations, given the 
% inferred statistic. Entropy levels are averaged across simulations and 
% saved.
%
% The results should be compared with Falk (1975)
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

addpath ../Tools/
addpath ../IdealObserversCode/

% Ideal observer options
in.jump         = 1;                % estimate with jumps
in.mode         = 'HMM';            % use the HMM (not sampling) algorithm
n               = 20;               % resolution of the univariate probability grid
in.verbose      = 1;                % check that no default values are used.


% Parameters of the simulation
pAlt = 0.05:0.05:0.95;
pJump = [1/3 1/6 1/12 1/24 1/96];  % a priori probability that a jump occurs at each outcome
L     = 21;                 % length of the sequence
nSim  = 1e4;			    % number of simulations

% Initialize variables
nAlt = numel(pAlt);
npJump = numel(pJump);
avg_entropy_trans = nan(nAlt, npJump, nSim);
avg_entropy_freq = nan(nAlt, npJump, nSim);
avg_entropy_pAlt = nan(nAlt, npJump, nSim);

% Run simulations
for iSim = 1:nSim
    fprintf('\n Comparison with jump, iter %d/%d', iSim, nSim)
    for ipAlt = 1:nAlt
        for ipJump = 1:npJump
            
            % generate a sequence sAlt of alternation and repetition,
            % with exactly the number of alt/rep specified.
            nAlt = round(pAlt(ipAlt)*(L-1));
            nRep = round((1-pAlt(ipAlt))*(L-1));
            sAlt = [1*ones(1, nAlt), 2*ones(1, nRep)];
            
            % randomize the order of rep/alt in this sequence
            sAlt = sAlt(randperm(L-1));
            
            % translate sAlt into a sequence of stimuli that alternate and
            % repeat
            s = zeros(1,L); s(1) = 1;
            for k = 2:L;
                if sAlt(k-1) == 1; s(k) = ~s(k-1); else s(k) = s(k-1); end
            end
            s = s+1;
            
            % compute posterior estimates given transition probabilities
            in.s            = s;                
            in.learned      = 'transition';     % estimate transition
            in.opt.pgrid    = linspace(0,1,n);  % univariate probability grid for estimation
            in.opt.Alpha0   = ones(n)/(n^2);    % uniform prior on transition probabilities
            in.opt.pJ       = pJump(ipJump);	% prior on volatility
            est_trans       = IdealObserver(in);
            
            % compute posterior estimates given stimulus frequency
            in.s            = s;                
            in.learned      = 'frequency';      % estimate transition
            in.opt.pgrid    = linspace(0,1,n);  % univariate probability grid for estimation
            in.opt.Alpha0   = ones(n,1)/n;      % uniform prior on transition probabilities
            in.opt.pJ       = pJump(ipJump);	% prior on volatility
            est_freq        = IdealObserver(in);
            
            % compute posterior estimates given probability to alternate or
            % repeat. To achieve this, the sequence is first converted into
            % "alternate or repeat", and then we use the observer that 
            % learns the frequency of outcomes (here, of alternations).
            in.s            = abs(diff(s))+1;
            in.learned      = 'frequency';      % estimate transition
            in.opt.pgrid    = linspace(0,1,n);  % estimation probability grid
            in.opt.Alpha0   = ones(n,1)/n;      % uniform prior on transition probabilities
            in.opt.pJ       = pJump(ipJump);
            est_pAlt        = IdealObserver(in);
            
            
            % compute entropy of these posterior estimates
            H_trans = - est_trans.p1_mean(end)*log2(est_trans.p1_mean(end)) ...
                - (1-est_trans.p1_mean(end))*log2(1-est_trans.p1_mean(end));
            H_freq = - est_freq.p1_mean(end)*log2(est_freq.p1_mean(end)) ...
                - (1-est_freq.p1_mean(end))*log2(1-est_freq.p1_mean(end));
            H_pAlt = - est_pAlt.p1_mean(end).*log2(est_pAlt.p1_mean(end)) ...
                - (1-est_pAlt.p1_mean(end)).*log2(1-est_pAlt.p1_mean(end));
            
            % compute average entropy
            avg_entropy_trans(ipAlt, ipJump, iSim) = mean(H_trans);
            avg_entropy_freq(ipAlt, ipJump, iSim) = mean(H_freq);
            avg_entropy_pAlt(ipAlt, ipJump, iSim) = mean(H_pAlt);
        end
    end
end

% average over simulations
m_avg_entropy_trans = mean(avg_entropy_trans, 3);
m_avg_entropy_freq = mean(avg_entropy_freq, 3);
m_avg_entropy_pAlt = mean(avg_entropy_pAlt, 3);

% save the results
timestamp = sprintf('%d-%d-%d_%d-%d-%1.0f_%5.0f', clock, rand*1e5);
save(['Comparison_NonStationarity_data_exactProb_',timestamp, '.mat'], ...
    'in', 'pAlt', 'pJump', 'nSim', ...
    'm_avg_entropy_trans', 'm_avg_entropy_freq', 'm_avg_entropy_pAlt')
