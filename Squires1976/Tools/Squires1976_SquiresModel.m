function [ p1_mean, surprise ] = Squires1976_SquiresModel( s, pA )
%SQUIRES1976_SQUIRESMODEL implements the linear model of Squires et al.
%(1976) published in Science. Note that they fitted it only on the last
%stimulus of each pattern whose length equals 5 (e.g. A in BBBBA).
%   - "s": the sequence of stimuli on which to run the model.
%   - "pA": the global probability (a scalar), that is magically provided.
%   - "Ndepth": the length of patterns on which to exert the local
%   probability pass (it is to 5 by Squires et al.).
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

%% INITIALIZATION
%  ==============

% Translate the sequence with 0s and 1s
if numel(unique(s)) > 2, error('The sequence is not binary'); end
s = [s, min(s)]; % add a A at the end of the sequence
news = NaN(size(s));
news(s == 1) = 1; % A
news(s == 2) = 0; % B
s = news;
N = numel(s);
Ndepth = 5;

% Fitted parameters
alpha = 0.6;
betaM = 0.235;
betaA = 0.033;
betaP = 0.505;

% Prepare output
p1_mean = NaN(size(s));

% For each stimulus
for k = Ndepth:N
    
    % Get the current chunk
    chunk = s(k-Ndepth+1:k-1);
    
    %% LOCAL PROBABILITY PASS
    %  ======================
    
    % Exponentially decaying memory
    forgetting = alpha .^ (Ndepth-1:-1:1);
    
    % We want the surprise evoked by the last stimulus
    M = sum(forgetting .* chunk);
    
    %% LOCAL ALTERNATION PASS
    %  ======================
    
    % Count the number of alternations
    if     ismember(chunk(end-3:end), [1 0 1 0], 'rows'), A = +3;
    elseif ismember(chunk(end-3:end), [0 1 0 1], 'rows'), A = -3;
    elseif ismember(chunk(end-3:end), [0 0 1 0], 'rows'), A = +2;
    elseif ismember(chunk(end-3:end), [1 1 0 1], 'rows'), A = -2;
    else A = 0;
    end
    
    %% COMPUTE THE EXPECTANCY SCORE
    %  ============================

    p1_mean(k) = (betaM * M) + ...  % local probability
                 (betaA * A) + ...  % local alternations
                 (betaP * pA) + ... % global probability
                 - 0.027;           % offset
end

% Compute surprise
p1_mean  = p1_mean(2:end);
surprise = ComputeSurprise(p1_mean, s(1:end-1));

% Compute Shannon surprise given predictions and actual outcomes.
function surp = ComputeSurprise(p1_mean, s)
seqL = numel(s);
surp = nan(1, seqL);
for j = 2:seqL
    if s(j) == 1
        surp(j) = -log2(p1_mean(j-1));      % likelihood of s(k)=1 given s(1:k-1)
    else
        surp(j) = -log2(1-p1_mean(j-1));    % likelihood of s(k)=2 given s(1:k-1)
    end
end
end

end
