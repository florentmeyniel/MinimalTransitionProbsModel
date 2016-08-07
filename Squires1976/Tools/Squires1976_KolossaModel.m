function [ p1_mean, surprise ] = Squires1976_KolossaModel( s )
%SQUIRES1976_KOLOSSAMODEL implements the model proposed by Kolossa et al.
%(2013) in Frontiers in Human Neuroscience.
%   - "s": the binary sequence
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

%% INTIALIZATION
%  =============

% Define useful variables %%% delete this??? k is not an argument of the
% function.....
% ~~~~~~~~~~~~~~~~~~~~~~~

% Get the counts of one of the stimuli
k = 1; % As
g = (s == k);

% The fitted free parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
alphaL  = 0.83;
tau1    = 33.6;
tau2    = 0.27;
alphaS  = 0.12;
betaS   = 1.82;
alphaD  = 0.05;
gammaD2 = 0.94;

% DECAY FACTOR FOR THE UPDATE OF LONG/SHORT TERM MEMORY
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The short-term memory is updated iteratively by forgetting past observartions.
% The decay factor is constant, corresponding to a memory that decays
% exponentially.
gammaS = exp(-1 / betaS);

% The long-term memory is udpate online. The current stimulus is
% integretaged into the long-term memory. The relative weight between a 
% given observation and the long-term estimate change over time, since the
% impact a one single observation become more and more marginal as more
% observations were accumulated. 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
betaL  = @(n) exp(+(((1/tau1) * n) + (1/tau2))); % NB: with a + sign here (typo in equation A2 in the paper??)
gammaL = @(n) exp(-1 / betaL(n));


% WEIGHTS OF PREVIOUS OBSERVATIONS FOR THE DETECTION OF ALTERNATIONS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Alternations are detected in patterns restricted to 4 observations.
% Observation are "convolved with an alternation kernel gammaD.
% Note that the alternated sign in the kernels detect values. Also note 
% that there is decay within this kernel: a recent alternation is weighted
% more that a remote alternation within the recent history
gammaDmax = 1;
gammaD = [gammaDmax-gammaD2, -(gammaDmax-gammaD2), gammaD2, -gammaD2]; % correspond to [gamma_4, gamma_3, gamma_2, gamma_1]
   

% NORNALIZING CONSTANTS
% ~~~~~~~~~~~~~~~~~~~~~
% These constants ensure that the result is in the range [0 1], i.e. that
% it is a probability.
% The additive normalizing contant
C  = 2; % correspond to (gammaDmax + gammaD(2) + gammaD(4)) / gammaDmax;

% The multiplicative normalizing constant
CD = 2*gammaDmax;

%% Run an iterative estimation of the generative probability
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% append a dummy stimulus at the end, so that the last observation of the 
% sequence provided as an input to this function is taken into account for 
% the inference of the hidden probability.
g = [g(:); 1]';

% initialize
p1_mean = NaN(1,numel(g));

% Starting values
% ~~~~~~~~~~~~~~~
cS_kn = 1/2; % unbiased prior for a binary sequence 
cL_kn = 1/2; % unbiased prior for a binary sequence

p1_mean(1) = (alphaL * cL_kn) + ...        % short-term memory
          (alphaS * cS_kn) + ...        % long-term memory
          (alphaD * (0 + (1/C)));       % alternation expectation

for n = 2:numel(g)
    
    % Update short term memory
    % ========================
    cS_kn = (1 - gammaS) * g(n-1) + (gammaS * cS_kn);
    
    % Update long term memory
    % =======================
    cL_kn = (1 - gammaL(n-1)) * g(n-1) + (gammaL(n-1) * cL_kn); 
    
    % Compute an alternation score given the recent history
    % =====================================================
    if n >= 5
        cD_kn = g(n-4:n-1)*gammaD' / CD;
    else
        cD_kn = NaN; % skip the first values
    end
    
    % COMBINE ALL THESE COMPONENTS INTO A PROBABILITY ESTIMATE
    % ========================================================
    p1_mean(n) = (alphaL * cL_kn) + ...        % short-term memory
              (alphaS * cS_kn) + ...        % long-term memory
              (alphaD * (cD_kn + (1/C)));   % alternation expectation
          
end

% In this code (like in Kolossa et al), p1_mean(t) is the probability
% of observing the item "1", given the previous observations s(1) to s(t-1)
% In our other codes, p1_mean(t) should be the probability of observing the
% item "1" given the previous observation s(1) to s(t). 
% Therefore, we shift the output for the sake of consistency with our other
% model.
p1_mean = p1_mean(2:end);
surprise = ComputeSurprise(p1_mean, s);

% Compute Shannon surprise given predictions and actual outcomes.
function surp = ComputeSurprise(p1_mean, s)
seqL = numel(s);
surp = nan(1, seqL);
surp(1) = -log2(0.5);
for k = 2:seqL
    if s(k) == 1
        surp(k) = -log2(p1_mean(k-1));      % likelihood of s(k)=1 given s(1:k-1)
    else
        surp(k) = -log2(1-p1_mean(k-1));    % likelihood of s(k)=2 given s(1:k-1)
    end
end
end

end
