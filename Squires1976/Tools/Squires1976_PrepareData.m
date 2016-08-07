function [ Data, Npat, p1list, PatLabels ] = Squires1976_PrepareData( DataFile )
%SQUIRES1976_PREPAREDATA loads Squires et al. (1976) or Kolossa et al.
%(2013) data and generates a matrix with the "tree-format".
%   - "DataFile": the name of the *.mat file that contains the data.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

%% FIND WHETHER WE SHOULD LOAD RESULTS FROM SQUIRES OR KOLOSSA
%  ===========================================================

% Get only the file name in case the whole path is specified
fs = strfind(DataFile, filesep);
if ~isempty(fs), FileName = DataFile(fs(end)+1:end);
else FileName = DataFile; end

% Guess which paper it is
if     any(strfind(FileName, 'Squires')), paper = 1;
elseif any(strfind(FileName, 'Kolossa')), paper = 2;
else error('The input data file is not recognized.');
end

%% GET DATA HEADER
%  ===============

% Load the data
OriginalData = load(DataFile);

% Get the global probabilities in force in each block
if     paper == 1, p1list = OriginalData.pA;
elseif paper == 2, p1list = OriginalData.pX;
end

% Get dimensions
Npat   = max(max(cellfun(@(x) numel(x), OriginalData.label)));
MaxLen = size(ff2n(Npat),1)/2;
Seq    = AllSeqPattern(Npat, 1);
Np1    = numel(p1list);

% If there is a single list, it means that it is the same across blocks
% (such as in Squires)
if any(size(OriginalData.label) == 1)
    OriginalData.label = repmat(OriginalData.label(:), [1,Np1]);
end

%% GET P300 AMPLITUDE
%  ==================

% For each p(A)
Data = NaN(Npat, MaxLen, Np1);
PatLabels = cell(Npat, MaxLen);
for iP = 1:Np1
    
    % Create a patterns' list sorted by size and by alphabetical order
    l = OriginalData.label(:,iP); % block's labels
    nl = {};
    for i = Npat:-1:1, nl = [nl; sort(l(cellfun(@numel, l) == i))]; end
    
    % Reorder the values according to that well-sorted list
    v = OriginalData.values(:,iP);
    [~,idx] = ismember(nl,l);
    l = l(idx);
    v = v(idx);
    
    % Discover whether patterns are based on A or B. In both case they will
    % be translater by associating a 1 to the letter on which they are
    % based (e.g. BABA: 2121 whereas AABB: 2211). The last letter
    % determines how the patterns are build and, by convention, we
    % attribute a 1 to that letter and a 2 to the other note.
    % in Kolossa et al. paper.
    seqlength = cellfun(@numel, nl);
    letter = nl(seqlength == 1);
    
    % Convert this new list from A & B to 1 & 2
    for i = 1:numel(seqlength)
        tmpseq = 2*ones(1,seqlength(i));
        tmpseq(strfind(l{i}, letter{1})) = 1;
        col = find(ismember(Seq{seqlength(i)}, tmpseq, 'rows'));
        Data(seqlength(i),col,iP) = v(i);
        PatLabels{seqlength(i),col} = l{i};
    end
end

end
