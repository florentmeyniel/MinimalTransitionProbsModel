function [MAP, m_h, s_h] = bern_ComputeMAPandPrediction(...
    NA, NB, pNA, pNB)
% Compute posterior MAP, mean and standard deviation, using the analytical formula for beta distributions.
%
% The advantage of using conjugate distribution is that posterior estimates
% (mean, variance, MAP) have analytical solutions that depend purely on the
% event counts NX augmented by the prior event counts pNX.
%
% Copyright 2016 Florent Meyniel & Maxime Maheu

% Beta parameters are equals to the event counts + 1.
% We convert event counts into beta parameters:
NA = NA + 1;
NB = NB + 1;

% The prior and the likelihood are both beta functions; their product is another
% beta function, whose parameters are the sum of paramaters of each beta
% distribution - 1
% We thus compute the parameters of the posterior beta distributions as:
NA = NA + pNA - 1;
NB = NB + pNB - 1;

% USE ANALYTICAL SOLUTION FOR BETA FUNCTIONS
MAP = (NA-1)/(NA+NB-2);
m_h = NA/(NA+NB);
s_h = sqrt(m_h * (1-m_h) / (NA + NB + 1));

end
