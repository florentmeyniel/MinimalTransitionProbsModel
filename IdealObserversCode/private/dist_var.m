function v = dist_var(val, pgrid, doWarning)
% Compute the variance (v) of a probability distribution (val) over a discrete 
% grid (pgrid). 
%
% Usage:
% v = dist_var(val, pgrid, doWarning)
% NB: doWarning: 1 allow a warning message on whether the distribution is 
% normalized to sum to 1 (default: 0) 
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
% check that a grid is provided
if nargin == 1
    error('pgrid is missing');
end
if nargin == 2
    doWarning = 0;
end

% Get number of distribution & distrubtion resolution
ndist = size(val, 1);
nres  = size(val, 2);

% check that val is a probability distribution
if ndist == 1
    if abs(sum(val) - 1) > 2*eps
        msg = sprintf(['the distribution is not a probability distribution! '...
            '\n ... it is now normalized.']);
        val = val / sum(val);
    end
    if any(abs(sum(val, 2) - 1) > 2*eps)
        msg = sprintf(['one or more distribution is not a probability distribution! '...
            '\n ... it is now normalized.']);
        val = val ./ repmat(sum(val, 2), [1, nres]);
    end
end

if doWarning == 1
    dist(msg)
end

% Get the mean
m = dist_mean(val, pgrid);

% compute variance
if ndist == 1
    v = sum(pgrid.^2 .* val) - m^2;
else
    if size(pgrid, 1) == 1
        % replicate grid for parallel estimation
        pgrid = repmat(pgrid, [ndist, 1]);
    end
        
    v = sum(pgrid.^2 .* val, 2) - m.^2;
end
