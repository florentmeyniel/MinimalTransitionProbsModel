function KL = dist_KLdiv(P, Q, doWarning)
% Compute the Kullback–Leibler divergence between probability distributions 
% P and Q. Note that KL is not symmetric. KL(Post, Prior) quantifies the
% information gain when moving from the prior to the posterior.
% 
% Usage:
% KL = dist_KLdiv(P, Q, doWarning)
% NB: doWarning: 1 allow a warning message on whether the distribution is 
% normalized to sum to 1 or the divergence not defined (default: 0) 
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
if nargin == 2
    doWarning = 0;
end

% Check that both probability distribution sum to 1
if abs((sum(Q) - 1)) > 2*eps
    msg = sprintf(['the 1st distribution is not a probability distribution! '...
            '\n ... it is now normalized.']);
    Q = Q / sum(Q);
end
if abs((sum(P) - 1)) > 2*eps
    msg = sprintf(['the 2nd distribution is not a probability distribution! '...
        '\n ... it is now normalized.']);
    P = P / sum(P);
end

% check absolute continuity (when Q is 0, then P is 0)
indQ0 = find(Q == 0);
indP0 = find(P == 0);
if length(intersect(indP0, indQ0)) ~= length(indQ0)
    msg = sprintf(['KL not defined because at least one value satifies:', ...
        'Q(i) = 0 & P(i) ~= 0']);
    KL = NaN;
else
    % When P(i) is zero, the contribution is interpreted as zero (since lim
    % of xlog(x) = 0 when x approaches 0). So these terms are simply
    % removed from the sum.
    ind = P > 0;
    KL = sum(P(ind) .* (log(P(ind)) - log(Q(ind))), 2);
end

if doWarning == 1
    dist(msg)
end
