% This script simulates many sequences with different apparent 
% probabilities of alternation. The ideal observer inference is computed 
% assuming dynamic belief. We distinguish between learning transition 
% probabilities, learning the stimulus frequency and learning the 
% frequency of alternations. We also compare different values 
% of the leak factor in the integration (from weak to strong). The entropy
% of the last posterior prediction is computed to quantify the 
% unpredictability of the sequence of observations, given the inferred 
% statistic. Entropy levels are averaged across simulations and saved.
%
% The results should be compared with Falk (1975)
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

addpath ../Tools/
addpath ../IdealObserversCode/

% Ideal observer options
in.jump             = 0;                % estimate assuming no jump
in.opt.AboutFirst   = 'WithoutFirst';   % discard 1st observation for analytical solutions when learning transition probabilities
in.opt.ReturnDist   = 0;                % Do not return full posterior distributions
in.opt.priorp1g2    = [1 1];            % uniform Beta prior
in.opt.priorp2g1    = [1 1];            % uniform Beta prior
in.opt.priorp1      = [1 1];            % uniform Beta prior
in.verbose          = 1;                % check that no default values are used.

% Parameters of the simulation
pAlt = 0.05:0.05:0.95;                  % generative probability of alternation
memDecay = [3 6 12 24 48 96];           % Memory decay
nSim     = 1e4;                         % number of simulations
L        = 21;

% Initialize variables
nAlt = numel(pAlt);
nD = numel(memDecay);
avg_entropy_trans = nan(nAlt, nD, nSim);
avg_entropy_freq = nan(nAlt, nD, nSim);
avg_entropy_pAlt = nan(nAlt, nD, nSim);

% Run simulations
for iSim = 1:nSim
    fprintf('\n Comparison without jump, iter %d/%d', iSim, nSim)
    for ipAlt = 1:nAlt
        for iD = 1:nD
            
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
            
            % Set memory of the ideal observer
            in.opt.MemParam     = {'Decay', memDecay(iD)};                  
            
            % compute posterior estimates, given transition probabilities
            in.s = s;
            in.learned = 'transition';
            est_trans = IdealObserver(in);
            
            % compute posterior estimates, given stimulus frequency
            in.s = s;
            in.learned = 'frequency';
            est_freq = IdealObserver(in);
            
            % compute posterior estimates, given probability of
            % alternations. To achieve this, the sequence is first
            % converted into "alternate or repeat", and then we use the
            % observer that learns the frequency of outcomes (here, of
            % alternations). Note that after conversion, the sequence is 1
            % element shorter. This does not bias the comparison since the
            % first element is not taken into account in the observer that
            % learns transition probabilities (cf. in.opt.AboutFirst)
            in.s = abs(diff(s))+1;
            in.learned = 'frequency';
            est_pAlt = IdealObserver(in);
            
            % compute entropy of these posterior estimates
            H_trans = - est_trans.p1_mean(end)*log2(est_trans.p1_mean(end)) ...
                - (1-est_trans.p1_mean(end))*log2(1-est_trans.p1_mean(end));
            H_freq = - est_freq.p1_mean(end)*log2(est_freq.p1_mean(end)) ...
                - (1-est_freq.p1_mean(end))*log2(1-est_freq.p1_mean(end));
            H_pAlt = - est_pAlt.p1_mean(end).*log2(est_pAlt.p1_mean(end)) ...
                - (1-est_pAlt.p1_mean(end)).*log2(1-est_pAlt.p1_mean(end));
            
            % compute average entropy
            avg_entropy_trans(ipAlt, iD, iSim) = mean(H_trans);
            avg_entropy_freq(ipAlt, iD, iSim) = mean(H_freq);
            avg_entropy_pAlt(ipAlt, iD, iSim) = mean(H_pAlt);
        end
    end
end

% average over simulations
m_avg_entropy_trans = mean(avg_entropy_trans, 3);
m_avg_entropy_freq = mean(avg_entropy_freq, 3);
m_avg_entropy_pAlt = mean(avg_entropy_pAlt, 3);

% save the results
timestamp = sprintf('%d-%d-%d_%d-%d-%1.0f_%5.0f', clock, rand*1e5);
save(['Comparison_Stationarity_data_exactProb_',timestamp, '.mat'], ...
    'in', 'pAlt', 'L', 'memDecay', 'nSim', ...
    'm_avg_entropy_trans', 'm_avg_entropy_freq', 'm_avg_entropy_pAlt')

