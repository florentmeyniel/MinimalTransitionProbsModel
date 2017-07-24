function out = IdealObserver(in)
% This function is a wrapper to run different Ideal Observer models on a
% sequence of binary values (with 1s and 2s). The models depend on:
%   WHAT IS LEARNED: frequency of the outcome, or transition probabilities 
%       between successive outcomes
%   WHETHER JUMPS ARE EXPECTED: the Ideal observer may assume that the
%       statistics that it estimates are stationary or change over time
%       with a given (fixed) probability at each new observation
%   HOW JUMPS ARE COMPUTED: In principle, two methods are available, one 
%       that uses a Hidden Markov chain Model (HMM) and one that uses 
%       sampling (see Meyniel et al Plos Comp Biol 2015). Here, we provide
%       only the code for HMM which is much faster.
%   WHETHER MEMORY IS PERFECT OR NOT in the case where jumps are not
%       expected. We provide two options that are not mutually exclusive:
%       memory is limited to a fixed sliding window, or assign weight to
%       past event with a exponential decay (= "leaky integration").
%
% Usage: out = IdealObserver(in)
%   
%   in: structure with the following fields:
%    in.learned: 'transition' or 'frequency'
%    in.jump: 1 if jumps are expected and 0 otherwise
%    in.mode: 'HMM' (default)
%    in.opt.pgrid: probability grid for numeric computations. It should be
%       a vector of values in the interval [0 1]. It is used when analytical
%       solutions are not available (HMM and algorithms without
%       jumps that do not assume that the 1st observation is arbitrary).
%       The default is 0:0.01:1
%    in.opt.WhichPass: 'Forward' (default and unique solution here, see the
%       HMM codes for more flexibility).
%    in.opt.priorp1: prior belief on the frequency of outcome 1, expressed as
%       the parameter of a Beta distribution. Default = [1 1]
%    in.opt.priorp1g2: prior belief on the transition probability p(1|2),
%       expressed as the parameter of a Beta distribution. Default [1 1].
%    in.opt.priorp2g1: prior belief on the transition probability p(2|1),
%       expressed as the parameter of a Beta distribution. Default [1 1].
%    in.opt.method: 'slow' (default) to compute posterior estimates given
%       each observation or 'quick' to compute only given all observations. 
%       This methods does not apply to HMM algorithms.
%    in.opt.AboutFirst: this parameter is used for the algorithm that expect
%       no jump and estimates transition probabilities. If 'WithoutFirst', 
%       (defaults) the algorithm assumes that the first outcome is 
%       arbitrary; in this case, the estimation is analytical. If 
%       'WithFirst', the 1st sample is assumed to depend on transition 
%       probabilities, which prevents the use of analytical solutions. 
%       Instead, numeric estimation of a grid is used. 
%    in.opt.MemParam: memory parameter when jumps are not expected.
%           [] (default) for a perfect memory
%           {'Limited', 10} a memory limited to the 10 most recent samples
%               (this number can be any given integer)
%           {'Decay', 2} a memory with decay exp(-n/2) applied to event n
%               in the past. (2 or any other integer).
%           {'Limited', 10, 'Decay', 2} for a combination.
%    in.opt.ReturnDist: this parameter is used for algorithms that expect
%       no jump. If set to 0, the full posterior distributions are not
%       return. They are if set to 1, but this requires some extra
%       computations.
%    in.verbose: 1 to point flag when default values are used
% 
%   out: structure with the following fields as results.
%    out.FunctionName: name of the Matlab function that implement the Ideal
%       Observer model estimation.
%    out.p1_mean: prediction (= expected likelihood) that 1 is the next outcome
%    out.p1_sd: estimation uncertainty (standard deviation) associated 
%       to this numeric prediction.
%    out.p1_dist: full posterior distribution when a frequency is estimated
%    out.p1g2_dist: full posterior distribution when p(1|2) is estimated
%    out.p2g1_dist: full posterior distribution when p(2|1) is estimated
%    out.p1g2xp2g1: full joint posterior distribution when the model
%       estimates p(1|2) and p(2|1)
%    out.surprise: Shannon surprise associated to the actual outcome, given the prediction.
%    out.Update: model update induced by each observation. This is computed
%       as a Kullback Leibler divergence on the posterior distributions.
% 
%    Note that all the output variables are time series. Element k in a
%    time series corresponds to the posterior estimate made the
%    observations 1 to k included.
% 
% Use IdealObserver('defaults') to return the defaults values.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
% ASSIGN DEFAULT VALUES
% =====================
indefaults.mode             = 'HMM';            % prefer HMM over sampling to estimate jumps
indefaults.opt.pgrid        = 0:0.01:1;         % estimation grid
indefaults.opt.WhichPass    = 'Forward';        % only choice here
indefaults.opt.priorp1g2    = [1 1];            % Bayes-Laplace prior
indefaults.opt.priorp2g1    = [1 1];            % Bayes-Laplace prior
indefaults.opt.priorp1      = [1 1];            % Bayes-Laplace prior
indefaults.opt.method       = 'slow';           % compute for each observation
indefaults.opt.MemParam     = [];               % perfect integration
indefaults.opt.AboutFirst   = 'WithoutFirst';   % to have analytical solutions
indefaults.opt.ReturnDist   = 0;                % do not return distribution when not already computed

% BASIS CHECKS ON INPUT
% =====================
if ischar(in)
    if strcmpi(in, 'defaults')
        out = indefaults;
        return
    else
        error('unknown input type')
    end
else
    in = CheckInput(in, indefaults, 'general');
end


% RUN CHECKS AND ESTIMATES FOR EACH TYPE OF OBSERVER
% ==================================================
if      strcmpi(in.learned, 'transition') && in.jump == 1 && strcmpi(in.mode, 'HMM')
    
    % ---------------------------------------------------------------------
    %         HMM ESTIMATION OF TRANSITION PROBABILITIES AND JUMPS
    % ---------------------------------------------------------------------
    
    out.FunctionName = 'ForwardBackward_MarkovJump';
    
    % CHECK INPUTS
    % ~~~~~~~~~~~~
    in = CheckInput(in, indefaults, 'pJ');
    in = CheckInput(in, indefaults, 'pgrid');
    in = CheckInput(in, indefaults, 'Alpha0 transition');
    in = CheckInput(in, indefaults, 'WhichPass');
    
    % RUN THE IDEAL OBSERVER
    % ~~~~~~~~~~~~~~~~~~~~~~
    % get the forward estimation of transition probabilities
    [~, rAlpha, ~, ~, Trans] = ForwardBackward_MarkovJump(in.s, in.opt.pJ, ...
        in.opt.pgrid, in.opt.Alpha0, in.opt.WhichPass);
    
    % translate the posterior distribution into a predictive distribution, 
    % which takes into account the assumed non-stationarity of the 
    % underlying statistics. 
    PredDist = TurnPosteriorIntoPrediction(rAlpha, Trans, in.opt.pJ);
        
    % Return full posterior distribution and marginal distributions
    out.p1g2xp2g1  = PredDist;
    out.p1g2_dist  = squeeze(sum(PredDist, 2));
    out.p2g1_dist  = squeeze(sum(PredDist, 1));
    
    % get "posterior" prediction (= mean) and confidence in this prediction
    % given each transition type
    out.p1g2_mean  = dist_mean(out.p1g2_dist', in.opt.pgrid)';
    out.p1g2_sd    = dist_sd(  out.p1g2_dist', in.opt.pgrid)';
    out.p2g1_mean  = dist_mean(out.p2g1_dist', in.opt.pgrid)';
    out.p2g1_sd    = dist_sd(  out.p2g1_dist', in.opt.pgrid)';
    
    % get prediction and confidence given the current observation
    out.p1_mean    = GetConditionalValue(out.p1g2_mean, 1-out.p2g1_mean, in.s);
    out.p1_sd      = GetConditionalValue(out.p1g2_sd,   out.p2g1_sd,     in.s);
    
    % get surprise given these predictions and the actual outcomes
    out.surprise   = ComputeSurprise(out.p1_mean, in.s);
    
    % get update of the posterior distribution (KL divergence)
    % NB: add the prior as 1st estimate.
    out.distUpdate  = DistributionUpdateMarkov(cat(3, in.opt.Alpha0, rAlpha));
    
    
elseif  strcmpi(in.learned, 'transition')  && in.jump == 1 && strcmpi(in.mode, 'sampling')
    
    % ---------------------------------------------------------------------
    %     GIBBS SAMPLING ESTIMATION OF TRANSITION PROBABILITIES AND JUMP
    % ---------------------------------------------------------------------
    out.FunctionName = 'Sampling_MarkovJump';
    
elseif  strcmpi(in.learned, 'transition') && in.jump == 0
    
    % ---------------------------------------------------------------------
    %         ESTIMATION OF TRANSITION PROBABILITIES WITHOUT JUMPS
    % ---------------------------------------------------------------------
    
    out.FunctionName = 'MarkovNoJump';
    
    % CHECK INPUTS
    % ~~~~~~~~~~~~
    in = CheckInput(in, indefaults, 'priorp1g2');
    in = CheckInput(in, indefaults, 'priorp2g1');
    in = CheckInput(in, indefaults, 'method');
    in = CheckInput(in, indefaults, 'MemParam');
    in = CheckInput(in, indefaults, 'AboutFirst');
    in = CheckInput(in, indefaults, 'ReturnDist');
    if in.opt.ReturnDist == 0 && strcmpi(in.opt.AboutFirst, 'WithoutFirst')
        % in this case, pgrid is not used
        in.opt.pgrid = [];
    else
        in = CheckInput(in, indefaults, 'pgrid');
    end
    
    
    % RUN THE IDEAL OBSERVER
    % ~~~~~~~~~~~~~~~~~~~~~~
    % get the forward estimation of transition probabilities
    [~, Marginal_mean, Marginal_sd] = MarkovNoJump(in.s, in.opt.pgrid, ...
        in.opt.priorp1g2, in.opt.priorp2g1, in.opt.method, ...
        in.opt.MemParam, [], in.opt.AboutFirst);
    
    % get "posterior" prediction (= mean) and confidence in this prediction
    % given each transition type
    out.p1g2_mean = Marginal_mean(1,:);
    out.p1g2_sd   = Marginal_sd(1,:);
    out.p2g1_mean = Marginal_mean(2,:);
    out.p2g1_sd   = Marginal_sd(2,:);
    
    % get prediction and confidence given the current observation
    out.p1_mean    = GetConditionalValue(out.p1g2_mean, 1-out.p2g1_mean, in.s);
    out.p1_sd      = GetConditionalValue(out.p1g2_sd,   out.p2g1_sd,     in.s);
    
    % get surprise given these predictions and the actual outcomes
    out.surprise   = ComputeSurprise(out.p1_mean, in.s);
    
    % Compute equivalent beta parameters
    [p1g2_a, p1g2_b] = BetaMomentsToParams(out.p1g2_mean, out.p1g2_sd);
    [p2g1_a, p2g1_b] = BetaMomentsToParams(out.p2g1_mean, out.p2g1_sd);
    
    % Return full posterior distribution and marginal distributions
    if in.opt.ReturnDist == 1
        out.p1g2_dist = zeros(numel(in.opt.pgrid), numel(in.s));
        out.p2g1_dist = zeros(numel(in.opt.pgrid), numel(in.s));
        out.p1g2xp2g1 = zeros(numel(in.opt.pgrid), numel(in.opt.pgrid), numel(in.s));
        for k = 1:numel(in.s)
            % beta distribution
            out.p1g2_dist(:,k) = betapdf(in.opt.pgrid, p1g2_a(k), p1g2_b(k));
            out.p2g1_dist(:,k) = betapdf(in.opt.pgrid, p2g1_a(k), p2g1_b(k));
            
            % normalize to have a discretized probability distribution
            out.p1g2_dist(:,k) = out.p1g2_dist(:,k) / sum(out.p1g2_dist(:,k));
            out.p2g1_dist(:,k) = out.p2g1_dist(:,k) / sum(out.p2g1_dist(:,k));
            
            % compute joint distribution
            out.p1g2xp2g1(:,:,k) = out.p1g2_dist(:,k) * out.p2g1_dist(:,k)';
        end
    end
    
    % get update of the posterior distribution (KL divergence)
    % (here, we use analytical solution for independent Beta distributions)
    % NB: add the prior as 1st estimate.
    out.distUpdate  = UpdateBetaJoint(...
        cat(2, in.opt.priorp1g2', [p1g2_a; p1g2_b]), ...
        cat(2, in.opt.priorp2g1', [p2g1_a; p2g1_b]));
    
elseif  strcmpi(in.learned, 'frequency')  && in.jump == 1 && strcmpi(in.mode, 'HMM')
    
    % ---------------------------------------------------------------------
    %                 HMM ESTIMATION OF FREQUENCIES AND JUMPS
    % ---------------------------------------------------------------------
    
    out.FunctionName = 'ForwardBackward_BernoulliJump';
    
    % CHECK INPUTS
    % ~~~~~~~~~~~~
    in = CheckInput(in, indefaults, 'pJ');
    in = CheckInput(in, indefaults, 'pgrid');
    in = CheckInput(in, indefaults, 'Alpha0 frequency');
    in = CheckInput(in, indefaults, 'WhichPass');
    
    % RUN THE IDEAL OBSERVER
    % ~~~~~~~~~~~~~~~~~~~~~~
    % get the forward estimation of frequencies
    [~, rAlpha, ~, ~, Trans] = ForwardBackward_BernoulliJump(in.s, in.opt.pJ, ...
        in.opt.pgrid, in.opt.Alpha0, in.opt.WhichPass);
    
    % translate the posterior distribution into a predictive distribution, 
    % which takes into account the assumed non-stationarity of the 
    % underlying statistics. 
    PredDist = TurnPosteriorIntoPrediction(rAlpha, Trans, in.opt.pJ);
    
    % Return full posterior distribution and marginal distributions
    out.p1_dist = PredDist;
    
    % get "posterior" prediction (= mean) and confidence in this prediction
    % given each transition type
    out.p1_mean = dist_mean(out.p1_dist', in.opt.pgrid)';
    out.p1_sd   = dist_sd(  out.p1_dist', in.opt.pgrid)';
    
    % get surprise given these predictions and the actual outcomes
    out.surprise = ComputeSurprise(out.p1_mean, in.s);
    
    % get update of the posterior distribution (KL divergence)
    % NB: add the prior as 1s estimate.
    out.distUpdate  = DistributionUpdateBernoulli(cat(2, in.opt.Alpha0, rAlpha));
    
    
elseif  strcmpi(in.learned, 'frequency')  && in.jump == 1 && strcmpi(in.mode, 'sampling')
    
    % ---------------------------------------------------------------------
    %       GIBBS SAMPLING ESTIMATION OF FREQUENCIES AND JUMP
    % ---------------------------------------------------------------------
    out.FunctionName = 'Sampling_BernoulliJump';
    
elseif  strcmpi(in.learned, 'frequency')  && in.jump == 0
    
    % ---------------------------------------------------------------------
    %         ESTIMATION OF FREQUENCIES WITHOUT JUMPS
    % ---------------------------------------------------------------------
    
    out.FunctionName = 'BernoulliNoJump';
    
    % CHECK INPUTS
    % ~~~~~~~~~~~~
    in = CheckInput(in, indefaults, 'priorp1');
    in = CheckInput(in, indefaults, 'method');
    in = CheckInput(in, indefaults, 'MemParam');
    in = CheckInput(in, indefaults, 'ReturnDist');
    if in.opt.ReturnDist == 1
        in = CheckInput(in, indefaults, 'pgrid');
    end
    
    
    % RUN THE IDEAL OBSERVER
    % ~~~~~~~~~~~~~~~~~~~~~~
    % get the forward estimation of frequencies (mean and confidence)
    [~, out.p1_mean, out.p1_sd] = BernoulliNoJump(in.s, in.opt.priorp1, ...
        in.opt.method, in.opt.MemParam, []);
    
    % get surprise given these predictions and the actual outcomes
    out.surprise = ComputeSurprise(out.p1_mean, in.s);
    
    % Compute equivalent beta parameters
    [p1_a, p1_b] = BetaMomentsToParams(out.p1_mean, out.p1_sd);
    
    % get update of the posterior distribution (KL divergence)
    % (here, we use analytical solution for Beta distributions)
    % NB: add the prior as 1s estimate.
    out.distUpdate  = UpdateBeta(cat(2, in.opt.priorp1', [p1_a; p1_b]));
    
    % Return full posterior distribution
    if in.opt.ReturnDist == 1
        out.p1_dist = zeros(numel(in.opt.pgrid), numel(in.s));
        for k = 1:numel(in.s)
            % beta distribution
            out.p1_dist(:,k) = betapdf(in.opt.pgrid, p1_a(k), p1_b(k));
            
            % normalize to have a discretized probability distribution
            out.p1_dist(:,k) = out.p1_dist(:,k) / sum(out.p1_dist(:,k));
        end
    end
end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % SUBFUNCTIONS % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

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

% Select the value conditional on the current observation
function val = GetConditionalValue(val_g2, val_g1, s)
seqL = numel(s);
val = nan(1, seqL);
for k = 1:seqL
    if s(k) == 1
        val(k) = val_g1(k);
    else
        val(k) = val_g2(k);
    end
end
end

% Return the KL divergence between two consecutive 2D distributions
function KLdiv = DistributionUpdateMarkov(Fdist)
% Initialize variables
seqL  = size(Fdist, 3);
KLdiv = nan(1, seqL);

% Convert 0 into eps (the closest value to 0 given Matlab largest
% double-precision) to compute the KL divergence without the numeric
% problem of 0*log(0), which occurs when absolute continuity does not hold.
% NB: These 0s (if any) are on the edges of the distribution. An
% alternative is to crop the distribution and discard the edge values (for
% a 2D distribution, this corresponds to discarding the first and last rows
% and the first and last columns). Both options, cropping or converting 0s
% to eps, asymptotically converge when the resolution of the grid % increases.
Fdist_corr = Fdist;
Fdist_corr(Fdist == 0) = eps; % approximate 0 value with Matlab largest double-precision

% Compute KL divergence across successive posterior estimates
for k = 2:seqL
    % Posterior distribution given stim 1 to k-1
    dist_prev = squeeze(Fdist_corr(:, :, k-1));
    
    % Posterior distribution given stim 1 to k
    dist_curr = squeeze(Fdist_corr(:, :, k));
    
    % KL divergence, i.e. distribution Update (note that the KL is not symmetric
    % so the order matters!! KL(Post || Prior) quantify how much
    % information is gained when one revises Prior into Post)
    KLdiv(k) = dist_KLdiv(dist_curr(:)', dist_prev(:)');
end

% Note that the first element in Fdist is actually the prior. Therefore, 
% currently, KLdiv(1) is nan, KLdiv(2) is the update between prior and 1st
% posterior, KLdiv(3) is the update between 1st and 2nd posterior, etc. We
% just want to get rid of the 1st nan and align KLdiv(k) on the update
% induced by the k-th observation.
KLdiv = KLdiv(2:end);
end

% Return the KL divergence between two consecutive 1D distributions
function KLdiv = DistributionUpdateBernoulli(Fdist)
% Initialize variables
seqL  = size(Fdist, 2);
KLdiv = nan(1, seqL);

% Convert 0 into eps (the closest value to 0 given Matlab largest
% double-precision) to compute the KL divergence without the numeric
% problem of 0*log(0), which occurs when absolute continuity does not hold.
% NB: These 0s (if any) are on the edges of the distribution (first and 
% value). Both options, cropping or converting 0s to eps, asymptotically 
% converge when the resolution of the grid increases.
Fdist_corr = Fdist;
Fdist_corr(Fdist == 0) = eps; % approximate 0 value with Matlab largest double-precision

% Compute KL divergence across successive posterior estimates
for k = 2:seqL
    % Posterior distribution given stim 1 to k-1
    dist_prev = squeeze(Fdist_corr(:, k-1));
    
    % Posterior distribution given stim 1 to k
    dist_curr = squeeze(Fdist_corr(:, k));
    
    % KL divergence, i.e. distribution Update (note that KL is not symmetric
    % so the order matters!! KL(Post || Prior) quantifies how much
    % information is gained when one revises Prior into Post
    KLdiv(k) = dist_KLdiv(dist_curr(:)', dist_prev(:)');
end

% Note that the first element in Fdist is actually the prior. Therefore, 
% currently, KLdiv(1) is nan, KLdiv(2) is the update between prior and 1st
% posterior, KLdiv(3) is the update between 1st and 2nd posterior, etc. We
% just want to get rid of the 1st nan and align KLdiv(k) on the update
% induced by the k-th observation.
KLdiv = KLdiv(2:end);
end

% Analytical KL divergence for Beta distributions
function KLdiv = UpdateBeta(all_param)
% Initialize variables
seqL  = size(all_param, 2);
KLdiv = nan(1, seqL);

% Compute KL divergence across successive posterior estimates
for k = 2:seqL
    % Get Beta parameters before and after the update
    param_post  = all_param(:,k)';
    param_prior = all_param(:,k-1)';
    
    % Use analytical solution for KL(post || prior) when post and prior
    % are Beta (or more generally Dirichlet) distributions.
    % cf e.g. http://bariskurt.com/kullback-leibler-divergence-between-two-dirichlet-and-beta-distributions
    % or https://en.wikipedia.org/wiki/Beta_distribution
    try
        KLdiv(k) = gammaln(sum(param_post)) - gammaln(sum(param_prior)) ...
            - sum(gammaln(param_post)) + sum(gammaln(param_prior)) ...
            + (param_post - param_prior) * (psi(param_post) - psi(sum(param_post)))';
    catch
        % for some special prior, the beta distribution is ill-defined.
        KLdiv(k) = nan;
    end
end

% Note that the first element in all_param is actually the prior. Therefore, 
% currently, KLdiv(1) is nan, KLdiv(2) is the update between prior and 1st
% posterior, KLdiv(3) is the update between 1st and 2nd posterior, etc. We
% just want to get rid of the 1st nan and align KLdiv(k) on the update
% induced by the k-th observation.
KLdiv = KLdiv(2:end);
end

% Analytical KL divergence for joint Beta distributions
function KLdiv = UpdateBetaJoint(all_param_d1, all_param_d2)
% Note: KL is additive for independent distributions.
KLdiv1 = UpdateBeta(all_param_d1);
KLdiv2 = UpdateBeta(all_param_d2);
KLdiv = KLdiv1 + KLdiv2;
end


% Equivalent parameters for a Beta distribution with mean m and sd s
function [a, b] = BetaMomentsToParams(m, s)
% Apply formula
v = s.^2;
a = ((1-m)./v -1./m).*(m.^2);
b = a.*(1./m - 1);

% round up to a given precision to avoid degeneracy
a = round(a*1e5)/1e5;
b = round(b*1e5)/1e5;
end

% Translate the posterior into a prediction by taking into account the
% assumed non stationarity
function outDist = TurnPosteriorIntoPrediction(inDist, TransMatrix, pJ)
% Initialize variables
seqL    = size(inDist);
seqL    = seqL(end);
outDist = nan(size(inDist));

% apply the assumed non-stationarity to the prediction
if numel(size(inDist)) == 3
    marg_length = size(inDist,1);
    tmpOut = nan(marg_length*marg_length, size(inDist,3));
    tmpIn  = reshape(inDist, marg_length*marg_length, size(inDist,3));
    for k = 1:seqL
        tmpOut(:,k) = ...
            (1-pJ)*tmpIn(:,k) + ...             % when there is no change
            pJ * (TransMatrix' * tmpIn(:,k));   % when there is a change
    end
    outDist = reshape(tmpOut, marg_length, marg_length, size(inDist,3));
elseif numel(size(inDist)) == 2
    for k = 1:seqL
        outDist(:,k) = ...
            (1-pJ)*inDist(:,k) + ...            % when there is no change
            pJ * (TransMatrix' * inDist(:,k));  % when there is a change
    end
else
    error('inDist should have 2 or 3 dimensions, not %d', inDist)
end
end

% Perform various checks on inputs
function in = CheckInput(in, indefaults, action)
switch action
    
    case 'general'
        if isstruct(in)
            if ~isfield(in, 'verbose')
                in.verbose = 1;
            else
                if ~ismember(in.verbose, [0 1])
                    error('unknown verbose mode')
                end
            end
            if ~isfield(in, 'learned')
                error('in should have a field ''learned''')
            end
            if ~isfield(in, 'jump')
                error('in should have a field ''jump''')
            end
            if ~ismember(lower(in.learned), lower({'transition', 'frequency'}))
                error('in.learned should be ''transition'' or ''frequency''')
            end
            if isnumeric(in.jump)
                if in.jump == 1
                    if ~isfield(in, 'mode')
                        in.mode = indefaults.mode;
                        if in.verbose; fprintf('\n use default for ''mode'''); end
                    end
                    if ~ismember(lower(in.mode), lower({'HMM', 'sampling'}))
                        error('in.mode should be ''HMM'' or ''sampling''')
                    end
                elseif in.jump ~= 0
                    error('in.jump should be equal to 0 or 1')
                end
            else
                error('in.jump should be equal to 0 or 1')
            end
        else
            error('unknown input type')
        end
        if ~isfield(in, 's')
            error('in should have a field ''s'' for the sequence of observations')
        else
            if ~isnumeric(in.s)
                error('in.s should be numeric')
            else
                if ~all(ismember(unique(in.s), [1 2]))
                    error('in.s should contain only 1s and 2s')
                else in.s = in.s(:)'; % force to a vector
                end
            end
        end
        if ~isfield(in, 'opt')
            error('in should have a field ''opt'' to specifies options')
        end
        
    case 'pJ'
        if ~isfield(in.opt, 'pJ')
            error(['a field in.opt.pJ is required for the HMM estimation of ', ...
                'transition probabilities with jumps'])
        else
            if ~isnumeric(in.opt.pJ) && (in.opt.pJ < 0 || in.opt.pJ > 1)
                error('in.opt.pJ should be a numeric value in the interval [0 1]')
            end
        end
        
    case 'pgrid'
        if ~isfield(in.opt, 'pgrid')
            in.opt.pgrid = indefaults.opt.pgrid;
            if in.verbose; fprintf('\n use default for ''pgrid'''); end
        else
            in.opt.pgrid = sort(unique(in.opt.pgrid));
            if min(in.opt.pgrid) < 0 || max(in.opt.pgrid) > 1
                error('in.opt.pgrid should be a vector of values in the interval [O 1]')
            end
        end
        
    case 'Alpha0 transition'
        if ~isfield(in.opt, 'Alpha0')
            in.opt.Alpha0 = ones(numel(in.opt.pgrid))/(numel(in.opt.pgrid))^2;
            if in.verbose; fprintf('\n use default for ''Alpha0'''); end
        end
        
    case 'Alpha0 frequency'
        if ~isfield(in.opt, 'Alpha0')
            in.opt.Alpha0 = ones(numel(in.opt.pgrid),1)/numel(in.opt.pgrid);
            if in.verbose; fprintf('\n use default for ''Alpha0'''); end
        else
            in.opt.Alpha0 = in.opt.Alpha0(:);
        end
        
    case 'WhichPass'
        if ~isfield(in.opt, 'WhichPass')
            in.opt.WhichPass = indefaults.opt.WhichPass;
        else
            if ~strcmpi(in.opt.WhichPass, 'Forward')
                error(['the only possible pass in.opt.WhichPass is ''Forward''', ...
                    'to run the Backward (and combined) pass, call the Ideal', ...
                    'Observer function directly, without this wrapper'])
            end
        end
        
    case 'priorp1g2'
        if ~isfield(in.opt, 'priorp1g2')
            in.opt.priorp1g2 = indefaults.opt.priorp1g2;
            if in.verbose; fprintf('\n use default for ''priorp1g2'''); end
        else
            if ~isnumeric(in.opt.priorp1g2) && numel(in.opt.priorp1g2) ~= 2
                error('in.opt.priorp1g2 should be a vector with two numbers');
            end
        end
        
    case 'priorp2g1'
        if ~isfield(in.opt, 'priorp2g1')
            in.opt.priorp2g1 = indefaults.opt.priorp2g1;
            if in.verbose; fprintf('\n use default for ''priorp2g1'''); end
        else
            if ~isnumeric(in.opt.priorp2g1) && numel(in.opt.priorp2g1) ~= 2
                error('in.opt.priorp2g1 should be a vector with two numbers');
            end
        end
        
    case 'method'
        if ~isfield(in.opt, 'method')
            in.opt.method = indefaults.opt.method;
        else
            if ~ischar(in.opt.method) && ~ismember(in.opt.method, {'slow', 'quick'})
                error('in.opt.method should be either ''slow'' or ''quick''')
            end
        end
        
    case 'MemParam'
        if ~isfield(in.opt, 'MemParam')
            in.opt.MemParam = indefaults.opt.MemParam;
            if in.verbose; fprintf('\n use default for ''MemParam'''); end
        else
            if ~(iscell(in.opt.MemParam) && any(length(in.opt.MemParam) == [2 4])) && ~isempty(in.opt.MemParam)
                error(['in.opt.MemParam should be either empty (perfect integration '...
                    'or a cell such as {''Limited'', 10} or  {''Decay'', 2}, or a combination ', ...
                    '{''Limited'', 10, ''Decay'', 2}'])
            elseif ~isempty(in.opt.MemParam)
                if length(in.opt.MemParam) == 2
                    if ~ismember(in.opt.MemParam{1}, {'Limited', 'Decay'})
                        error('in.opt.MemParam = %s is not recognized. It should be ''Limited'' or ''Decay''', ...
                            in.opt.MemParam{1})
                    end
                else
                    for k = 1:2
                        ind = [1 3];
                        if ~ismember(in.opt.MemParam{ind(k)}, {'Limited', 'Decay'})
                            error('in.opt.MemParam = %s is not recognized. It should be ''Limited'' or ''Decay''', ...
                                in.opt.MemParam{ind(k)})
                        end
                    end
                end
            end
        end
        
    case 'AboutFirst'
        if ~isfield(in.opt, 'AboutFirst')
            in.opt.AboutFirst = indefaults.opt.AboutFirst;
            if in.verbose; fprintf('\n use default for ''AboutFirst'''); end
        else
            if ~ischar(in.opt.AboutFirst) && ~ismember(in.opt.AboutFirst, {'WithoutFirst', 'WithFirst'})
                error('in.opt.AboutFirst should be either ''WithoutFirst'' or ''WithFirst''')
            end
        end
        
    case 'ReturnDist'
        if ~isfield(in.opt, 'ReturnDist')
            in.opt.ReturnDist = indefaults.opt.ReturnDist;
            if in.verbose; fprintf('\n use default for ''ReturnDist'''); end
        end
        
    case 'priorp1'
        if ~isfield(in.opt, 'priorp1')
            in.opt.priorp1 = indefaults.opt.priorp1;
            if in.verbose; fprintf('\n use default for ''priorp1'''); end
        else
            if ~isnumeric(in.opt.priorp1) && numel(in.opt.priorp1) ~= 2
                error('in.opt.priorp1 should be a vector with two numbers');
            end
        end
end
end
