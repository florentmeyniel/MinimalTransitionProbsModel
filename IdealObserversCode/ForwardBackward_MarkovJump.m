function [rGamma, rAlpha, rBeta, JumpPost, Trans] = ForwardBackward_MarkovJump(s, pJ, pgrid, Alpha0, Pass)
% Forward-Backward algorithm to solve a hidden markov model, in which the
% hidden states are non-stationary transition probabilities controlling
% the observed outcome s.
%
% Given that s is a sequence of binary outcomes (A or B), the hidden states
% are pairs of possible states (x1, x2) corresponding to independent
% transition probabilities (from A to B and from B to A).
%
% Usage:
% [rGamma, rAlpha, rBeta, JumpPost] = ...
%   ForwardBackward_MarkovJump(s, pJ, pgrid, Alpha0, WhichPass)
%
% Input:
%      s: sequence of numeric values (coded as 1s and 2s for convenience)
%     pJ: prior on jump occurrence (change in transition probabilities) at
%         any given moment (this is a scalar value)
%  pgrid: grid for numeric estimation of transition probabilities.
%         pgrid(i) is p(s(i)=1 | x1(i), s(i-1)=1)
%         The grid is the same for both transition probabilities.
% Alpha0: prior on states, i.e. on transition probabilities (pAgB varies
%         over rows & pBgA varies over columns)
%   Pass: if 'Forward' (default), computes only the forward estimation. If
%         'Backward', also computes the backward estimation and the
%         combination of forward and backward estimates.
%
% Output:
%   rGamma: rGamma(:,t) is the posterior distribution for states x(t)
%           (i.e. hidden transition probabilities) given observation s(1:N)
% 	    => the Forward-Backward estimation
%   rAlpha: rAlpha(:,t) is the posterior distribution for states x(t)
%           (i.e. hidden transition probabilities) given observation s(1:t)
% 	    => the Forward estimation
%    rBeta: rBeta(:,t) is the posterior distribution for states x(t)
%           (i.e. hidden transition probabilities) given observation s(t+1:N)
% 	    => the Backward estimation
% JumpPost: Posterior on jump probability given s(1:N).
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
% Checks and initialization
% =========================
seqL   = length(s);
n      = length(pgrid);
pgrid  = sort(pgrid(:))';   % make sure that pgrid is a row vector in ascending order
Alpha0 = Alpha0(:);         % make sure that Alpha0 is a column vector

if pJ > 1 || pJ < 0
    error('pJ must be in the [0 1] interval (current value: %d)', pJ)
end
if ~exist('Pass', 'var')
    Pass = 'Forward';
end

% FORWARD PASS: p(x1(i,t),x2(i,t) | s(1:t))
% =========================================

% Alpha(i,t) = p(x1(i,t), x2(i,t) | s(1:t))
% Since we want to distinguish jump vs. no jump, we split the values of
% Alpha in 2 columns, corresponding to jump or no jump.
% Alpha(i, 1, t) = p(x1(i,t), x2(i,t), J=0 | s(1:t))
% Alpha(i, 2, t) = p(x1(i,t), x2(i,t), J=1 | s(1:t))

% Initialize Alpha (forward estimation)
% Pairs of (x1, x2) are sorted into a column
% Transition (i,j) corresponds to: from x1(i) & x2(i) to x1(j) & x2(j)
Alpha = zeros(n*n, 2, seqL);
Alpha(:,1,1) = Alpha0/2;
Alpha(:,2,1) = Alpha0/2;

% get the matrix of non diagonal elements
NonDiag_nn = ones(n*n);
NonDiag_nn(logical(eye(n*n))) = 0;

% Compute the transition matrix (jump = non diagonal transitions).
% NB: the prior on jump occurrence pJ is NOT included here.
% Trans(i,j) is the probability to jump FROM i TO J
% Hence, sum(Trans(i,:)) = 1
% The likelihood of new values after a jump correspond the prior on states.
Trans = (1 ./ (NonDiag_nn * Alpha0)) * Alpha0';
Trans = Trans .* NonDiag_nn;

% Compute alpha iteratively (forward pass)
for t = 2:seqL;
    
    % Specify likelihood of current observation
    % LL(i, j) = p(s(t) | pAgB(i), pBgA(j))
    LL = zeros(n, n);
    if      s(t-1) == 1 && s(t) == 2
        LL = repmat(pgrid, [n, 1]);     % B|A
    elseif  s(t-1) == 1 && s(t) == 1
        LL = repmat(1-pgrid, [n, 1]);   % A|A
    elseif  s(t-1) == 2 && s(t) == 1
        LL = repmat(pgrid', [1, n]);    % A|B
    elseif  s(t-1) == 2 && s(t) == 2
        LL = repmat(1-pgrid', [1, n]);  % B|B
    end
    
    % Normalize LL so that it is a probability over the possible transition
    % values considered.
    LL = LL / sum(LL(:));
    
    % Sort LL as a column
    LL = LL(:);
    
    % Compute the new alpha, based on the former alpha, the prior on
    % transition between states and the likelihood. See for instance
    % Wikipedia entry for 'Forward algorithm'.
    if t == 2
        Alpha(:,1, t) = (1-pJ)   * LL .* Alpha0;
        Alpha(:,2, t) = pJ/(n-1) * LL .* Alpha0;
    else
        
        % No Jump at t:
        % - take the prior on 'no jump': (1-pJ)
        % - take the current observation likelihood under x_i (LL)
        % - take the posterior on x_i(t-1) (summed over whether there was a
        % jump of not at t-1)
        Alpha(:,1,t) = (1-pJ) * LL .* (...
            (Alpha(:,1,t-1) + Alpha(:,2,t-1)));
        
        % Jump at t:
        % - take the prior on 'jump': pJ
        % - take the current observation likelihood under x_i (LL)
        % - take the posterior on all the other states, excluding x_i(t-1)
        % (summed over whether there was a jump or not at i-1)
        % - sum over the likelihood of the ORIGINS of such state
        % (hence the transpose on Trans, to sum over the ORIGIN)
        Alpha(:,2,t) = pJ * LL .* (...
            Trans' * (Alpha(:,1,t-1) + Alpha(:,2,t-1)));
    end
    
    % scale alpha as a posterior (which we will do eventually) to avoid
    % numeric overflow
    NormalizationCst = sum(sum(Alpha(:,:,t),2),1);
    Alpha(:,1, t) = Alpha(:,1,t) / NormalizationCst;
    Alpha(:,2, t) = Alpha(:,2,t) / NormalizationCst;
end

% BACKWARD PASS: p(s(t+1:N | x(i,t), x(i,t), s(t))
% ================================================
if strcmpi(Pass, 'Backward')
    
    % Beta(i,t) = p(s(t+1:N) | x1(i,t), x2(i,t), s(t))
    % Since we want to distinguish jump vs. no jump, we split
    % the values of Beta in 2 columns, corresponding to jump or no jump.
    % Beta(i,1,t) = p(s(t+1:N), J=0 | x1(t), x2(t), s(t))
    % Beta(i,2,t) = p(s(t+1:N), J=1 | x1(t), x2(t), s(t))
    
    % In addition, we normalize Beta(i,t) so that it sums to 1 over i. This is
    % only for convenience (to make the interpretation of numeric values easier)
    % since in the end the backward probabilities are normalized, we can
    % apply any scaling factor to Beta(i,t)
    
    % Initialize beta (backward estimation)
    % Pairs of (x1, x2) are sorted into a column
    % Transition (i,j) corresponds to: from x1(i) & x2(i) to x1(j) & x2(j)
    Beta = zeros(n*n, 2, seqL);
    
    % Specify likelihood of pA and pB, depending on the transition
    % probabilities.
    % pA & pB are matrices corresponding to pAgB (change of rows) X pBgA (change of columns)
    pAgB = repmat(pgrid', [1 n]);
    pBgA = repmat(pgrid, [n 1]);
    
    pA = pAgB ./ (pAgB+pBgA);
    pB = 1-pA;
    
    % force these two to be 0.5 (not nan) to allow further computations
    pA(1,1) = 0.5;
    pB(1,1) = 0.5;
    
    % Compute beta iteratively (backward pass)
    for t = seqL:-1:2;
        
        if t == seqL
            Beta(:,1,t) = 1;
            Beta(:,2,t) = 1;
        else
            % Specify likelihood p(s(t)|x(i,t),s(i,t-1))
            LL = zeros(n, n);
            if      s(t-1) == 1 && s(t) == 2
                LL = repmat(pgrid, [n, 1]);     % B|A
            elseif  s(t-1) == 1 && s(t) == 1
                LL = repmat(1-pgrid, [n, 1]);   % A|A
            elseif  s(t-1) == 2 && s(t) == 1
                LL = repmat(pgrid', [1, n]);    % A|B
            elseif  s(t-1) == 2 && s(t) == 2
                LL = repmat(1-pgrid', [1, n]);  % B|B
            end
            
            % No Jump from t to t+1
            % Average over potential states x(t)
            % take only diagonal elements
            Beta(:,1, t) = (1-pJ) * LL(:) .* (...
                (Beta(:,1,t+1) + Beta(:,2,t+1)));
            
            % Jump from t to t+1
            % Average over potential states x(t)
            % sum over non diagonal elements
            % NB: there is no transpose here of Trans because we sum over
            % TARGET location (not ORIGIN)
            Beta(:,2, t) = pJ * (...
                 Trans * (LL(:) .*(Beta(:,1,t+1) + Beta(:,2,t+1))));
        end
        
        % scale beta to sum = 1. 
        NormalizationCst = nansum(nansum(Beta(:,:,t),2),1);
        Beta(:,1, t) = Beta(:,1,t) / NormalizationCst;
        Beta(:,2, t) = Beta(:,2,t) / NormalizationCst;
    end
    
    % Shift Beta so that Beta(:,:,t) is the posterior given s(t+1:N)
    newBeta = zeros(size(Beta));
    newBeta(:,:,1) = 1;
    newBeta(:,:,2:end) = Beta(:,:,1:end-1);
    clear Beta; Beta = newBeta; clear newBeta;
end

% COMBINE FORWARD AND BACKWARD PASS
% =================================
if strcmpi(Pass, 'Backward')
    % Compute the forward & backward posterior
    sBeta  = squeeze(sum(Beta, 2));
    rBeta  = reshape(sBeta, [n, n, seqL]);
    sAlpha = squeeze(sum(Alpha, 2));
    rAlpha = reshape(sAlpha, [n, n, seqL]);
    
    % p(x(t)|s(1:N)) ~ p(s(t+1:N)|x(t),s(t)) p(x(t)|s(1:t))
    % NB: the sum over the second dimension is to average out the J=0 or 1
    sGamma = sAlpha .*  sBeta;
    
    % Scale gamma as a posterior over observations
    cst = repmat(sum(sGamma, 1), [n*n, 1]);
    sGamma = sGamma ./ cst;
    rGamma = reshape(sGamma, [n, n, seqL]);
        
    % Compute the posterior on jump, summed over the states
    % GammaJ = p(x(t),J=1|s(1:N)) ~ p(x(t),J=1|s(1:t)) p(s(t+1:N)|x(t),s(t),J=1) 
    %                             ~ p(x(t),J=1|s(1:t)) [1/p(J=1) * p(s(t+1:N),J=1|x(t))] 
    %                             ~ Alpha2 [1/p(J=1) * Beta2] 
    GammaJ   = squeeze(Alpha(:,2,:)) .* ((1/pJ) * squeeze(Beta(:,2,:)));
    GammaJ = GammaJ ./ cst;
    JumpPost = sum(GammaJ, 1);
else
    rGamma = []; rBeta = []; JumpPost = [];
    sAlpha = squeeze(sum(Alpha, 2));
    rAlpha = reshape(sAlpha, [n, n, seqL]);
end

end
