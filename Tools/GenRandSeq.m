function [ s, gen_p1, gen_p1g2, gen_p2g1 ] = GenRandSeq( L, p )
%GENRANDSEQ is a simple function generating a random sequence from
%specified theoretical (events or transitions) probabilistic rules.
%   - "L": a Nx1 vector (with N being the number of chunks) specifying the
%   length of the chunks 
%   - "p": a Nx1 (Bernoulli) specifying p(1) in each chunk or a Nx2
%   (Markov) specifying p(1|2) and p(2|1) in each chunk.
%Usage:
%   >> s = GenRandSeq(100, 1/3); % will generate a sequence of 100 stimuli
%   with p(1) = 1/3
%   >> s = GenRandSeq(500, [1/3, 2/3]); % will generate a sequence of 500
%   stimuli with p(1|2) = 1/3 and p(2|1) = 2/3
%   >> s = GenRandSeq(100, 1/3); % will generate a sequence of 100 stimuli
%   with p(1) = 1/3
%   >> s = GenRandSeq([150, 100, 50], [1/3, 2/3, 1/3; 1/2, 1/2, 1/2]); will
%   generate a sequence composed of a 150 stimuli chunk with p(1|2) = 1/3 
%   and p(2|1) = 2/3, followed by a 100 stimuli chunk with p(1|2) = 1/2 and
%   p(2|1) = 1/2, followed by a 50 stimuli chunk with p(1|2) = 1/3 and
%   p(2|1) = 1/2
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% Check that the number of chunk length matches the number of
% probabilities' values
if numel(L) ~= size(p,1)
    error(['The number of chunks to make does not match ', ...
        'the number of probabilities'' values']);
end

% Initialize random generators
rng('shuffle')

% Number of chunks
Nchunks = numel(L);

% Chunks' limits
L = L(:)';
climits = [[1, cumsum(L(1:end-1))+1]', cumsum(L)'];

% Prepare sequence
s = NaN(1, sum(L));
gen_p1 = NaN(1, sum(L));
gen_p1g2 = NaN(1, sum(L));
gen_p2g1 = NaN(1, sum(L));

% For each chunk
for c = 1:Nchunks
    
    % Length of the chunk
    N = L(c);
    
    % Indices for the current chunk
    idx = climits(c,1):1:climits(c,2);
    
    % Events' probabilities (Bernoulli)
    if numel(p(c,:)) == 1
        
        % Theoretical probability
        p1 = p(c);
        gen_p1(idx)   = p1;
        gen_p1g2(idx) = p1;
        gen_p2g1(idx) = 1 - p1;
        
        % Generate the random sequence
        seq = rand(1, N);
        seq = seq(randperm(N));
        RndGen = seq - p1;
        s(idx) = (RndGen > 0) + 1;
        
    % Transitions' probabilities (Markov)
    elseif numel(p(c,:)) == 2
        
        % Theoretical probabilities
        p1g2 = p(c,1);
        p2g1 = p(c,2);
        p1   = p1g2 / (p1g2 + p2g1);
        gen_p1(idx)   = p1;
        gen_p1g2(idx) = p1g2;
        gen_p2g1(idx) = p2g1;
        
        % Transition matrix from [1 2] to 1
        T = [1 - p2g1, p1g2];
        
        % Generate the random sequence
        for k = 1:N
            if k == 1, RndGen = 1 + rand - p1;
            else       RndGen = rand - T(s(idx(k)-1)); end
            if     RndGen <= 0, s(idx(k)) = 1;
            elseif RndGen >  0, s(idx(k)) = 2; end
        end
    end
end

end
