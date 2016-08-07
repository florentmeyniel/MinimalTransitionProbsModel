function seq = AllSeqPattern( N, stim )
%ALLSEQPATTERN generates a list of all the possible patterns of size 1:N
%that end with a particular stimulus.
%   - "N": the maximum length of patterns.
%   - "stim": the stimulus (1,2,...) ending the patterns that have to be
%     generated.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% Default stimulus
if nargin < 2, stim = 1; end

% Find the patterns
seq = cell(1, N);
for k = 1:N
    seq{k} = [ff2n(k-1)+1, repmat(stim,2^(k-1),1)];
end

end
