function m = dist_mean(val, pgrid, doWarning)
% Compute the mean (m) of a probability distribution (val) over a discrete 
% grid (pgrid). 
%
% Usage:
% m = dist_mean(val, pgrid, doWarning)
%
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

if ndist == 1
    m = sum(val .* pgrid);
else
    % replicate grid for parallel estimation
    pgrid = repmat(pgrid, [ndist, 1]);
    m = sum(val .* pgrid, 2);
end
