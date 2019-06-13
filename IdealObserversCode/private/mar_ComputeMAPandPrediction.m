function [MAP, m_h, s_h] = mar_ComputeMAPandPrediction(NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB)
% Compute posterior MAP, mean and standard deviation, using the analytical formula for beta distributions.
%
% The advantage of using conjugate distribution is that posterior estimates
% (mean, variance, MAP) have analytical solutions that depend purely on the
% event counts NXgY (X given Y) augmented by the prior event counts pNygY.
%
% Copyright 2016 Florent Meyniel & Maxime Maheu

% Beta parameters are equals to the event counts + 1.
% We convert event counts into beta parameters:
NAgA = NAgA + 1;
NAgB = NAgB + 1;
NBgA = NBgA + 1;
NBgB = NBgB + 1;

% The prior and the likelihood are both beta functions; their product is another
% beta function, whose parameters are the sum of paramaters of each beta
% distribution - 1
% We thus compute the parameters of the posterior beta distributions as:
NAgA = NAgA + pNAgA - 1;
NAgB = NAgB + pNAgB - 1;
NBgA = NBgA + pNBgA - 1;
NBgB = NBgB + pNBgB - 1;

% USE ANALYTICAL SOLUTION
MAP = [(NAgB-1)/(NAgB+NBgB-2); (NBgA-1)/(NBgA+NAgA-2)];
m_h = [NAgB/(NAgB+NBgB); NBgA/(NBgA+NAgA)];

s_h = sqrt([m_h(1) * (1-m_h(1)) / (NAgB + NBgB + 1);...
    m_h(2) * (1-m_h(2)) / (NBgA + NAgA + 1)]);

end
