function TreeMat = ReformatChoDataAsSquiresTree(data, labels)
% This function reformats the data by Cho et al (2002) to plot them as a 
% tree like Squires et al (1976).
% 
% Usage: TreeMat = ReformatChoDataAsSquiresTree(data, labels)
%
% TreeMat is prepared so that it can then be plotted with PlotTree(TreeMat)
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% INITIALIZATION
% ===============

% convert patterns from R/A (repeating / alternation) to X/Y (one stimulus 
% or the other)
pattern = {
    'RRRR', 'XXXXX'; 
    'ARRR', 'YXXXX';
    'RARR', 'YYXXX';
    'AARR', 'XYXXX';
    'RRAR', 'YYYXX';
    'ARAR', 'XYYXX'; 
    'RAAR', 'XXYXX'; 
    'AAAR', 'YXYXX';
    'RRRA', 'YYYYX';
    'ARRA', 'XYYYX';
    'RARA', 'XXYYX';
    'AARA', 'YXYYX';
    'RRAA', 'XXXYX'; 
    'ARAA', 'YXXYX';
    'RAAA', 'YYXYX';
    'AAAA', 'XYXYX' 
    };

% Initialize output vector
TreeMat = zeros(5, 16);

% table to convert the binary output of ff2n into X/Y
Pat2Str = ['X', 'Y'];

% FILL DATA FROM CHO ET AL IN THE LAST COLUMN OF THE TREE
% =======================================================

% compute all possible combinations
pat = ff2n(4);

% loop over conditions
for iPat = 1:size(pat, 1);
    
    % translate pattern into string
    str = [];
    for k = 1:4;
        str = [str, Pat2Str(pat(iPat, k)+1)];
    end
    str = [str, 'X'];
    
    % get corresponding index in Cho et al
    ind = find(ismember({pattern{:,2}}, str));
    ind = find(ismember(labels, pattern{ind,1}));
    TreeMat(5, iPat) = data(ind);
end

% COMPUTE THE REMAINING NODES OF THE TREE WITH SUCCESSIVE AVERAGING
% =================================================================

for iL = 4:-1:1
    
    % compute number of node in this layer of the tree
    nPat = 2^(iL-1);
    
    % Average the patterns that differ only by the 1st element.
    % By construction with ff2n, these patterns are nPat rows apart in
    % TreeMat.
    for iPat = 1:nPat
        TreeMat(iL,iPat) = 0.5*(...
            TreeMat(iL+1, iPat) + ...
            TreeMat(iL+1, nPat+iPat)...
            );
    end
end
