function [varargout] = BernoulliNoJump(s, priorpA, method, MemParam, OutPut)
% Ideal Observer that estimates the probability (frequency) generating
% the observed outcomes.
%
% Usage:
% [MAP, ...
%   m_hat, ...
%   s_hat, ...
%   pmY, ...
%       ] = BernoulliNoJump(s, priorpA, method, MemParam, OutPut)
%
% Input:
%   * s is the sequence (1: A; 2: B)
%   * priorpA is the prior on p(A), expressed as the distribution of a beta
%     distribution: A ~ beta(priorpA(1), priorpA(2))
%   * method: 'slow' (default) to compute for every outcome, 'quick' to
%   compute only on the last.
%   * param: cell with paired arguments (name, value)
%       o the default is to use all events without decay
%       o {'Limited', 10}: the memory window size is limited to 10 events
%       o {'Decay', 2} an exponential decay exp(-1/2*n) is applied
%         to event n in the past. Note: the last even has n=1
%       o can be combined, e.g. {'Limited', 10, 'Decay, 2}
%   * OutPut: the output the only what is specified by OupPut. OutPut is a
%       string (output variable name) or a cell of strings for multiple
%       outputs.
%
% Output:
%   All output variables are time series (except in the 'quick' mode).
%   The k-th element in each output variable corresponds to the inferred
%   value given observations 1 to k s(1:k).
%
%   * MAP: maximum a posteriori probability
%   * m_hat: mean posterior probability
%   * s_hat: standard deviation of the posterior estimate
%   * pmY: log model evidence (marginal likelihood of observations received so far)
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
% PARSE INPUT
% ===========
% if no prior, use the Bayes-Laplace prior
if nargin == 2; priorpA = [1 1]; end

% A conjugate prior distribution is used for pA. This distribution is a
% Beta distributions; the parameters of a Beta distribution are equivalent
% to an event count -1.
% Here, we convert the (Beta) prior into an event count.
pNA = priorpA(1)-1;
pNB = priorpA(2)-1;

% Choose computation mode (time series, or given the entire sequence).
if nargin <= 2
    method = 'slow';
elseif ~(strcmp(method, 'quick') || strcmp(method, 'slow'))
    error('check method parameters')
end

% Set the memory parameter
if nargin <= 3
    MemParam = [];
elseif iscell(MemParam) && any(length(MemParam) == [2 4])
elseif isempty(MemParam)
else
    error('check memory parameters')
end

% Control output parameters
if nargin <= 4
    OutPut = [];
elseif ~exist('OutPut', 'var')
    error('check parameter indicating output')
end

% number of elements in the sequence
L = length(s);

% POSTERIOR INFERENCE
% ===================
if strcmp(method, 'slow')
    
    % MAKE AN ESTIMATION FOR EACH OBSERVATION IN THE SEQUENCE OF
    % OBSERVATION
    
    % Initialize variables
    MAP            = zeros(1, L);
    m_hat          = zeros(1, L);
    s_hat          = zeros(1, L);
    pmY            = zeros(L, 1);
    
    % FOR EACH OBSERVATIONS
    for k = 1:L
        
        % Compute event count, including the current observation
        [~, ~, ~, ~, ~, NA, NB] = CountEventInMemory(s(1:k), MemParam);
        
        % Compute MAP & predictions
        [MAP(k), m_hat(k), s_hat(k)] = bern_ComputeMAPandPrediction(NA, NB, pNA, pNB);
        
        
        % Compute Log Model Evidence
        % Since we use a conjugate prior, LL*Prior is a beta distribution,
        % and its integral can be computed analytically
        pmY(k) = log(beta(NA+pNA+1, NB+pNB+1));
        
        % Note that this is numerically equivalent (provided that the grid as
        % many data point) to:
        % Prior = (t.^pNA) .* ((1-t).^pNB);
        % LL = (t.^NA) .* ((1-t).^NB);
        % pmY = log((1/(nt-1))*sum(LL.*Prior))
        % Caution!! the 'dt' term in the integral should not be forgotten in
        % the discrete approximation! t here has 2 dimensions.
    end
    
else
    
    % MAKE AN ESTIMATION GIVEN ALL THE OBSERVED EVENTS
    % Compute event count, including the current observation
    [~, ~, ~, ~, ~, NA, NB] = CountEventInMemory(s, MemParam);
    
    % Compute MAP & predictions
    [MAP, m_hat, s_hat] = bern_ComputeMAPandPrediction(NA, NB, pNA, pNB);
    
    
    % Compute Log Model Evidence
    % Since we use a conjugate prior, LL*Prior is a beta distribution,
    % and its integral analytically
    pmY = log(beta(NA+pNA+1, NB+pNB+1));
end

% Control the output of the function
if isempty(OutPut)
    if nargout >= 1, varargout{1} = MAP;   end
    if nargout >= 2, varargout{2} = m_hat; end
    if nargout >= 3, varargout{3} = s_hat; end
    if nargout >= 4, varargout{4} = pmY;   end
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
