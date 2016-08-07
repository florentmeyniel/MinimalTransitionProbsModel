function [varargout] = MarkovNoJump(s, t, priorpAgB, priorpBgA, method, MemParam, OutPut, AboutFirst)
% Ideal Observer that estimates the hidden transition probabilities 
% generating the observed outomes. 
%
% Usage:
% [MAP, ...
%     m_hat, ... 
%     s_hat, ... 
%     pmY, ...   
%     predA, ...      
%     predA_sd, ...   
%           ] = MarkovNoJump(s, t, priorpAgB, priorpBgA, method, MemParam, OutPut, AboutFirst)
% Input:
%   * s is the sequence (1: A; 2: B)
%   * t is the (unidimensional) grid for numeric parameter estimation.
%   * priorpAgB is the prior on p(A|B), expressed as the distribution of a
%     beta distribution: A|B ~ beta(priorpAgB(1), priorpAgB(2))
%   * method: 'slow' (default) to compute on every outcome, 'quick' to
%   compute only on the last.
%   * param: cell with paired arguments (name, value)
%       o the default is to use all events without decay
%       o {'Limited', 10}: the memory window size is limited to 10 events
%       o {'Decay', 2} an exponential decay exp(-1/2*n) is applied 
%         to event n in the past. Note: the last even has n=1            % 
%       o can be combined, e.g. {'Limited', 10, 'Decay, 2}
%   * OutPut: the output the only what is specified by OupPut. OutPut is a 
%       string (output variable name) or a cell of strings for multiple 
%       outputs.
%   * AboutFirst: 'WithoutFirst', to assume (unlike the default option) 
%     that the 1st event is drawn randomly. In that case, analytic
%     posterior values can be computed.
%
% Ouput: 
%   All output variables are time series (except in the 'quick' mode). 
%   The k-th element in each output variable corresponds to the inferred 
%   value given observations 1 to k s(1:k). When the output variable has 2 
%   elements (MAP, m_hat, s_hat) they corresponds respectively to A|B and B|A.
% 
%   * MAP: maximum a posteriori transition probabilities, 
%   * m_hat: mean posterior transition probability
%   * s_hat: standard deviation of the posterior estimate
%   * pmY: log model evidence (marginal likelihood of observations received so far)
%   * predA: likelihood that the next observation is A, i.e. that s(k+1)=1
%   * predA_sd: standard deviation of the estimated likelihood of the next event
%
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 

% PARSE INPUT
% ===========
% if no prior, use the Bayes-Laplace prior
if nargin <= 3; priorpAgB = [1 1]; priorpBgA = [1 1]; end

% Conjugate prior distributions are used for pAgB and pBgA. These 
% distributions are Beta distributions; the parameters of a Beta 
% distribution are equivalent to an event count -1.
% Here, we convert the (Beta) prior into a event count.
pNAgB = priorpAgB(1)-1;
pNBgB = priorpAgB(2)-1;
pNBgA = priorpBgA(1)-1;
pNAgA = priorpBgA(2)-1;

% Choose computation mode (time series, or given the entire sequence).
if nargin <= 4
    method = 'slow';
elseif ~(strcmp(method, 'quick') || strcmp(method, 'slow'))
    error('check method parameters')
end

% Set the memory parameter
MemDecay = Inf; % events are counted with perfect, unlimited memory
if nargin <= 5
    MemParam = [];
elseif iscell(MemParam) && any(length(MemParam) == [2 4])
elseif isempty(MemParam)
else
    error('check memory parameters')
end

% Control output parameters
if nargin <= 6 
    OutPut = []; 
elseif ~exist('OutPut', 'var')
    error('check parameter indicating output')
end

% Specify the assumption about the first event
if nargin <= 7
    AboutFirst = 'WithFirst';
elseif ~(ischar(AboutFirst) && (strcmpi(AboutFirst, 'WithoutFirst') || strcmpi(AboutFirst, 'WithFirst')))
    error('Check ''AboutFirst'' argument')
end

% number of elements in the grid and the sequence
nt = length(t);
L  = length(s);

% POSTERIOR INFERENCE
% ===================
if strcmp(method, 'slow')
    
    % MAKE AN ESTIMATION FOR EACH OBSERVATION IN THE SEQUENCE OF
    % OBSERVATION
    
    % Initialize variables
    MAP            = zeros(2, L);
    m_hat          = zeros(2, L);
    s_hat          = zeros(2, L);
    predA          = zeros(1,L);
    predA_sd       = zeros(1,L);
    pmY            = zeros(L, 1);
    
    % FIRST OBSERVATION
    % The posterior estimate on transition probabilities after having 
    % observed the first event must be made before having received any 
    % observation about any transition: it therefore corresponds to the 
    % prior. 
    [MAP(:, 1), m_hat(:, 1), s_hat(:, 1)] = ComputeMAPandPrediction([], [], [], 'WithoutFirst', ...
        0, 0, 0, 0, pNBgA, pNAgB, pNAgA, pNBgB);
    
    % prediction about the next observation
    if s(1) == 1
        predA(1)    = 1-m_hat(2,1); % = 1-p(B|A)
        predA_sd(1) = s_hat(2,1);   
    else
        predA(1)    = m_hat(1,1); % = p(A|B)
        predA_sd(1) = s_hat(1,1);   
    end
    
    % SUBSEQUENT OBSERVATIONS
    for k = 2:L
        
        % Count transitions including the current one
        [NBgA, NAgB, NAgA, NBgB, s1] = CountEventInMemory(s(1:k), MemParam);
        
        % Likelihood of the 1st event
        if strcmpi(AboutFirst, 'WithoutFirst')
            % Assume that the first event followed an equiprobabable random
            % draw. This will afford analytical solutions.
            p1 = 1/2;
        else
            % Assume that the first event was sampled given the transition
            % probabilities p(A|B) and p(B|A). This will prevent analytical
            % solutions.
            p1 = ComputeLikelihood1stEvent(MemDecay, s1, k, nt, t);
        end
        
        % Compute the full posterior distribution only when the analytical 
        % solution is not available, i.e. only when the first event is
        % assumed to be governed by transitin probabilities 
        % (instead of being sampled equiprobabibly). 
        if ~strcmpi(AboutFirst, 'WithoutFirst')
            post = ComputePosterior(p1, t, NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB);
        else
            post = [];
        end
        
        % Compute MAP & predictions
        [MAP(:,k), m_hat(:,k), s_hat(:,k), predA(k), predA_sd(k)] ...
            = ComputeMAPandPrediction(post, s(1:k), t, AboutFirst, NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB);
        
        % Compute Log Model Evidence      
        if nargout >=4 || any(strcmp(OutPut, 'pmY')) % compute only if necessary because it is time consuming
            pmY(k) = ComputeLME(p1, nt, t, AboutFirst, NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB);
        end
        
    end
    
else
    
    % MAKE AN ESTIMATION GIVEN ALL THE OBSERVED EVENTS
    
    % Count transitions including the current one
    [NBgA, NAgB, NAgA, NBgB, s1] = CountEventInMemory(s, MemParam);
    
    % Likelihood of the 1st event
    if strcmp(AboutFirst, 'WithoutFirst')
        % Assume that the first event followed an equiprobabable random
        % draw. This will afford analytical solutions.
        p1 = 1/2;
    else
        % Assume that the first event was sampled given the transition
        % probabilities p(A|B) and p(B|A). This will prevent analytical
        % solutions.
        p1 = ComputeLikelihood1stEvent(MemDecay, s1, L, nt, t);
    end
    
    % Compute the full posterior distribution only when the analytical
    % solution is not available, i.e. only when the first event is
    % assumed to be governed by transitin probabilities
    % (instead of being sampled equiprobabibly).
    if ~strcmpi(AboutFirst, 'WithoutFirst')
        post = ComputePosterior(p1, t, NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB);
    else
        post = [];
    end
    
    % Compute MAP & predictions
    [MAP, m_hat, s_hat, predA, predA_sd] ...
        = ComputeMAPandPrediction(post, s, t, AboutFirst, NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB);
    
    % Compute Log Model Evidence
    if nargout >=4 || any(strcmp(OutPut, 'pmY')) % compute only if necessary because it is time consuming
        pmY = ComputeLME(p1, nt, t, AboutFirst, NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB);
    end
end

% Control the output of the function
if isempty(OutPut)
    if nargout >= 1, varargout{1} = MAP;      end
    if nargout >= 2, varargout{2} = m_hat;    end
    if nargout >= 3, varargout{3} = s_hat;    end
    if nargout >= 4, varargout{4} = pmY;      end
    if nargout >= 5, varargout{5} = predA;    end
    if nargout >= 6, varargout{6} = predA_sd; end
else
    if iscell(OutPut)
        for k = 1:length(OutPut)
            eval(sprintf('varargout{%d} = %s;', k, OutPut{k}))
        end
    else
        eval(sprintf('varargout{%d} = %s;', 1, OutPut))
    end
end

end

% % % % % % % % % % % % % % % % % SUBFUNCTIONS % % % % % % % % % % % % % %

function p1 = ComputeLikelihood1stEvent(MemDecay, s1, k, nt, t)
% 1: Compute likelihood of 1st event
% ----------------------------------
% This likelihood is expressed as a function of p(A|B) (columns) and p(B|A)
% (rows). In other words, it is expressed as a p(B|A) x p(A|B) matrix.
pA = (ones(nt, 1)*t) ./ (ones(nt, 1)*t + t'*ones(1,nt));
pB = (t'*ones(1, nt)) ./ (ones(nt, 1)*t + t'*ones(1,nt));

if s1 == 1 % first is A
    p1 = pA;
else
    p1 = pB;
end

% apply a decay on p1 that is the same as for the other (more recent)
% events. Note that if MemDecay = Inf, there is actually no decay.
p1 = p1.^(exp(-(1/MemDecay)*(k-1))); % -1: e.g. k=2: we have in memory events 1 & 2; 1 is 1 time back in the past.

end

function [post] = ComputePosterior(p1, t, NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB)
% COMPUTE POSTERIOR
% =================
% notations:
%   - t the transition proba
%   - y the data
%
% p(t|y) = p(y|t)*p(t)/p(y)
% p(y) does not depend on t.
% so p(t|y) ~ p(y|t)*p(t)
%
% p(y|t) = p(y1|t)
%           * [(1-t(0|1))^N(1|1)]*[t(0|1)^N(0|1)] ...
%           * [t(1|0)^N(1|0)]*[(1-t(1|0))^N(0|0)]
%
% note that t(0|1)|y ~ Beta(N(0|1)+1, N(1|1)+1)
%       and t(1|0)|y ~ Beta(N(1|0)+1, N(0|0)+1)


% 2: compute posterior for remaining events: p(y2,...,yn|t)*p(t)
% --------------------------------------------------------------
% use the beta distribution. Note that they are not normalized, which
% does not matter since only the shape of the distribution matters.
tmp1 = (betapdf(t, NBgA+pNBgA+1, NAgA+pNAgA+1))';
tmp2 = betapdf(t, NAgB+pNAgB+1, NBgB+pNBgB+1);

% 3: compute posterior
% --------------------
% post is pBgA (change of line) x pAgB (change of columns)
% In other words, it is expressed as a p(B|A) x p(A|B) matrix.
post = p1 .* (tmp1 * tmp2);

% normalize so that the histogram is a discrete probability
post = post ./ nansum(post(:), 1);

end

function [MAP, m_h, s_h, predA, predA_sd] = ComputeMAPandPrediction(post, s, t, AboutFirst, NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB)

if strcmpi(AboutFirst, 'WithoutFirst')
    % USE ANALYTICAL SOLUTION
    [MAP, m_h, s_h] = mar_ComputeMAPandPrediction(...
    NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB);

else    % USE ESTIMATION ON A GRID
    
    % 4: compute MAP & variance
    % ~~~~~~~~~~~~~~~~~~~~~~~~~
    [~, indpAgB] = max(max(post, [], 1));
    [~, indpBgA] = max(max(post, [], 2));
    
    MAP = t([indpAgB indpBgA]);
    
    % compute the marginal probabilities
    mpAgB = nansum(post, 2);
    mpBgA = nansum(post, 1)';
    
    m_h = [nansum(mpBgA.*t'); ...
        nansum(mpAgB.*t')];
    
    s_h = sqrt([nansum(t'.^2.*mpBgA) - m_h(1)^2; ...
        nansum(t'.^2.*mpAgB) - m_h(2)^2]);
end

% 5: compute likelihood of next event
% -----------------------------------
% p(y(n)|y(1,...,n-1)) = int[post * t(y(n), y(n-1)) dt]
% which is the posterior marginal expectancy of the transition rate
% computed just above.
if nargout > 3
    if s(end-1) == 1 % use posterior marginal exectancy of B|A
        predA = 1-m_h(2);
        predA_sd = s_h(2);
    else           % use posterior marginal exectancy of A|B
        predA = m_h(1);
        predA_sd = s_h(1);
    end
end
end

function pmY = ComputeLME(p1, nt, t, AboutFirst, NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB)
% Compute model evidence (the marginal likelihood)
% ------------------------------------------------
% p(y|Markov) = int[p(y|Markov,t)*p(t) dt]
% Note that p(y|Markov,t)*p(t) was computed above, but up to a scaling
% factor. Here, it is this scaling factor that matters.

% Values from the continuous pdf are computed, and integrated with a
% discrete approximation.
% Caution!! the 'dt' term in the integral should not be forgotten in
% the discrete approximation! t here has 2 dimensions.

if strcmpi(AboutFirst, 'WithoutFirst')
    % USE ANALYTICAL SOLUTION
    
    t1 = [NAgB + pNAgB + 1, NBgB + pNBgB + 1];
    t2 = [NBgA + pNBgA + 1, NAgA + pNAgA + 1];
    
    % log integral of a beta distribution with these parameters
    intBeta1 = sum(gammaln(t1)) - gammaln(sum(t1));
    intBeta2 = sum(gammaln(t2)) - gammaln(sum(t2));
    
    % The likelihood is the product of the two integrals
    LL = intBeta1 + intBeta2;
    
    % add the LL of the first event    
    pmY = log(p1) + LL;
    
else
    % USE ESTIMATION ON A GRID
    
    tmp1 = (betapdf(t, NAgA+1, NBgA+1)*beta(NAgA+1, NBgA+1))';
    tmp2 = (betapdf(t, NBgB+1, NAgB+1)*beta(NBgB+1, NAgB+1));
    LL = (p1 .* (tmp1 * tmp2));
    
    Prior = (betapdf(t, pNAgA+1, pNBgA+1)*beta(pNAgA+1, pNBgA+1))' * ...
        (betapdf(t, pNBgB+1, pNAgB+1)*beta(pNBgB+1, pNAgB+1));
    
    pmY = (1/(nt-1))^2*(nansum(nansum(LL .* Prior)));
    
    % change to log scale
    pmY = log(pmY);
    
end
end
