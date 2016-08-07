function [MAP, m_h, s_h] = bern_ComputeMAPandPrediction(...
    NA, NB, pNA, pNB)
% Compute posterior MAP, mean and standard deviation, using the analytical formula for beta distributions.
% 
% The advantage of using conjugate distribution is that posterior estimates
% (mean, variance, MAP) have analytical solutions that depend purely on the
% event counts NX augmented by the prior event counts pNX.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
% USE ANALYTICAL SOLUTION
MAP = (NA+pNA)/(NA+pNA+NB+pNB);
    
m_h = (NA+pNA+1)/(NA+pNA+NB+pNB+2);

s_h = sqrt( (NA+pNA+1)*(NB+pNB+1) / ((NA+pNA+NB+pNB+2)^2*(NA+pNA+NB+pNB+3)) );...

end
