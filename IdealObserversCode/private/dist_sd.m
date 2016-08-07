function s = dist_sd(val, pgrid)
% Compute the standard deviation (s) of a probability distribution (val) 
% over a discrete grid (pgrid). 
%
% Usage:
% s = dist_var(val, pgrid)
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
% check that a grid is provided
if nargin == 1
    error('pgrid is missing');
end

s = sqrt(dist_var(val, pgrid));
