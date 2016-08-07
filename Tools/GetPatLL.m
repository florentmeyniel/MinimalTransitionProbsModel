function [ mLL, nLL, sLL, seq, getLL, getLLpos, getLLcat, seqnames ] = ...
    GetPatLL( N, s, seqLL, sym, stim, encaps )
%GETPATLL seeks for and count patterns in a sequence besides
%returning the average "seqLL" quantity associated with each pattern.
%   - "N", the length of patterns to look for.
%   - "s", the sequence in which patterns have to be seek for.
%   - "seqLL", the variable to average for each patterns (it can be
%   theoretical surprise, event-related field value, ...).
%   - "sym", if it is true, then AAB is considered as BBA, the pattern
%   thus becomes "XXY".
%   - "stim" specifies the identity of the stimulus on which to focus on
%   (if the "sym" option is enable, this variable is not used).
%   - "encaps" specifies whether there should be a strict equivalence
%   between the number of patterns (i.e. A = AA+BA = AAA+ABA+BBA = ...).
%       * 0: seek for all patterns no matter of whether shorter ones should
%       be encapsulated in longer ones (default).
%       * 1: check only the NaN in the sequence (equivalent to "2" if there
%       is no NaN value in the "seqLL" vector).
%       * 2: check the NaN both in the sequence and in the variable to
%       average ("seqLL") to ensure completely encapsulated patterns.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

%% INITIALIZATION
%  ==============

% By default, seek for all of the possible patterns 
if nargin < 6
    encaps = false;
    
    % By default, we look for the A stimulus
    if nargin < 5
        stim = 1;
        
        % By default, A different from B
        if nargin < 4
            sym = false;
        end
    end
end

% Compute all sequence of stimuli
seq = AllSeqPattern(N);
L = numel(s);

% Initialize likelihood structure
getLL    = cell(size(seq, 2), size(seq{end}, 1));
getLLpos = cell(size(seq, 2), size(seq{end}, 1));
getLLcat = cell(size(seq, 2), size(seq{end}, 1));

% We do not care about which stimulus to look for if we consider XXY like
% patterns (A and B have a symmetrical role)
if sym == true, stim = 1; end

%% GET ALL THE PATTERNS
%  ====================
%  i.e. if there is  ..., NaN, 2, 1, ..., then nA = +1, nBA = +1 but no
%  3 stimuli-based patterns will be updated. Therefore there is no strict
%  dependancies between the number of observations in each pattern length.
if encaps == false

    % Loop over the events' numbers
    for kk = 1:N

        % Loop over the sequence to get the type of sequence the chunk of sequence belongs to
        for k = kk:L
            type = [];

            % Get the pattern present in this little sequence snip
            [type, cat] = GetPatternType(seq{kk}, s(k-kk+1:k), sym, stim);
            
            % If the type is recognized (i.e. if it does not end with 2 (B) in case of stim = 1)
            if ~isempty(type)
                getLL   {kk, type} = [getLL{kk, type}, seqLL(k)];
                getLLpos{kk, type} = [getLLpos{kk, type}, k];
                getLLcat{kk, type} = [getLLcat{kk, type}, cat];
            end
        end
    end

%% GET ONLY ENCAPSULATED PATTERNS
%  ==============================
%  i.e. nA = nAA + nBA = nAAA + nBBA + nABA ...
elseif encaps > 0
    
    % For each pattern length
    for kk = 1:N
        
        % For each stimulus from the pattern length to the sequence length
        for k = N:L
            type = [];
            
            % Get the current sequence
            seqpart = s(k-N+1:k);
            modpart = seqLL(k-N+1:k);
            
            % Depending on the number of conditions, the condition can be
            % more or less severe, the second one ensure completely
            % encapsulated patterns.
            if     encaps == 1, cond = sum(isnan(seqpart)) == 0;
            elseif encaps == 2, cond = sum(isnan(seqpart)) == 0 && sum(isnan(modpart)) == 0;
            end
            
            % If the max-length sequence if availble, we get the type of
            % pattern (it ensures that there is strict depencies between
            % patterns with different length)
            if cond
                [type, cat] = GetPatternType(seq{kk}, seqpart(end-kk+1:end), sym, stim);
            
                % If there was one of the patterns we are looking for in the
                % sequence snips, add the related value we are interested in
                if ~isempty(type)
                    getLL   {kk, type} = [getLL{kk, type}, seqLL(k)];
                    getLLpos{kk, type} = [getLLpos{kk, type}, k];
                    getLLcat{kk, type} = [getLLcat{kk, type}, cat];
                end
            end 
        end
    end
end

%% AVERAGE LIKELIHOOD FOR EACH PATTERN
%  ===================================

% Compute statistics
maxlen = size(ff2n(N),1)/2;
nLL = NaN(N,maxlen);
mLL = NaN(N,maxlen);
sLL = NaN(N,maxlen);
for kk = 1:N
    for k = 1:2^(kk-1)
        n = nansum(~isnan(getLL{kk,k}), 2); % Number of patterns
        if      isempty(n), nLL(kk,k) = 0;
        elseif ~isempty(n), nLL(kk,k) = n; end
        mLL(kk,k) = nanmean(getLL{kk,k}); % Mean likelihood of the pattern
        sLL(kk,k) = nanstd(getLL{kk,k}); % Standard-deviation likelihood of the pattern
    end
end

% Get sequences' names
if     sym == false
    if     stim == 1, stimnames = {'A', 'B'};
    elseif stim == 2, stimnames = {'B', 'A'};
    end
elseif sym == true,   stimnames = {'X', 'Y'};
end
seqnames = cell(size(seq));
for kk = 1:N
    seqnames{kk} = stimnames(seq{kk});
end

%% NESTED FUNCTION
%  ===============

function [ pattype, patcat ] = GetPatternType( patternslist, seqpart, reverse, stim )
% This function nicely seeks for any pattern specifices in "patternslist"
% that could be present in "seqpart". The other arguments are explained in
% the help section of the main function.

    % We try both patterns (e.g. AA is equal to BB)
    if reverse % we start by trying A-based patterns
        pattype = find(ismember(patternslist, seqpart, 'rows'));
        patcat  = 1;
        if isempty(pattype) % if no one is found, we try B-based patterns
            pattype = find(ismember(-patternslist+3, seqpart, 'rows'));
            patcat  = 2;
        end

    % We only try a single pattern (e.g. AA is different from BB)
    elseif ~reverse
        if     stim == 1 % A
            pattype = find(ismember( patternslist,   seqpart, 'rows'));
        elseif stim == 2 % B
            pattype = find(ismember(-patternslist+3, seqpart, 'rows'));
        end
        patcat = stim;
    end
end

end
