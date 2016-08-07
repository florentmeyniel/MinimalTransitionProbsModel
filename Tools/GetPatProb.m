function PatProb = GetPatProb( nPat, pA )
%PATPROB generates the theoretical probability of observing particular
%patterns given a certain p(A).
%   - "nPat": a scalar specifying the maximum length of patterns to
%   investigate.
%   - "pA": a scalar spcifying the p(A).
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% Get all the possible patterns
seq = AllSeqPattern(nPat);

% Prepare output
PatProb = NaN(nPat, size(ff2n(nPat),1)/2);

% For each pattern length
for k = 1:nPat
    tmp = zeros(size(seq{k}));

    % Get event expected probability
    tmp(seq{k} == 1) =   pA;
    tmp(seq{k} == 2) = 1-pA;

    % Compute the pattern expected probability
    w = prod(tmp, 2);
    PatProb(k,1:numel(w)) = w;
end

end
