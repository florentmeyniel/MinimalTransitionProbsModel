function [MAP, m_h, s_h] = mar_ComputeMAPandPrediction(NBgA, NAgB, NAgA, NBgB, pNBgA, pNAgB, pNAgA, pNBgB)
% Compute posterior MAP, mean and standard deviation, using the analytical formula for beta distributions.
% 
% The advantage of using conjugate distribution is that posterior estimates
% (mean, variance, MAP) have analytical solutions that depend purely on the
% event counts NXgY (X given Y) augmented by the prior event counts pNygY.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
% USE ANALYTICAL SOLUTION
MAP = [(NAgB+pNAgB)/(NAgB+pNAgB+NBgB+pNBgB); ...
    (NBgA+pNBgA)/(NBgA+pNBgA+NAgA+pNAgA)];
m_h = [(NAgB+pNAgB+1)/(NAgB+pNAgB+NBgB+pNBgB+2); ...
    (NBgA+pNBgA+1)/(NBgA+pNBgA+NAgA+pNAgA+2)];

s_h = sqrt([(NAgB+pNAgB+1)*(NBgB+pNBgB+1) / ((NAgB+pNAgB+NBgB+pNBgB+2)^2*(NAgB+pNAgB+NBgB+pNBgB+3));...
    (NBgA+pNBgA+1)*(NAgA+pNAgA+1) / ((NBgA+pNBgA+NAgA+pNAgA+2)^2*(NBgA+pNBgA+NAgA+pNAgA+3))]);

end
