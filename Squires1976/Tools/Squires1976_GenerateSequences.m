function Seq = Squires1976_GenerateSequences( pA, L, nSeq, nPat )
%SQUIRES1976_GENERATESEQUENCES generates "nSeq" sequences composed of "L"
%stimuli for each specified probability in the "pA" list, all composed with
%every possible patterns composed of "nPat" stimuli (i.e. if nPat = 4:
%AAAA, AAAB, ..., BBBA, BBBB).
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
% Get dimensions
nP = numel(pA);
Seq = NaN(nSeq, L, nP);

% Variables used for printing purpose
nd = numel(sprintf('%i', nSeq));
rp = 11+(2*nd);
bk = repmat('\b', [1,rp]);

% For each block and each sequence example
fprintf('\nGenerating the stimuli sequences...\n');
for blk = 1:nP
    fprintf(' - Block %i/%i [p(A) = %1.2f]:\t%s', blk, nP, pA(blk), repmat(' ', [1,rp]));
    for seq = 1:nSeq

        % Print the state
        fprintf([bk, 'sequence %', num2str(nd), '.0f/%' , num2str(nd), '.0f\n'], seq, nSeq);

        % While all patterns are not there, keep trying
        nMissing = +inf;
        while nMissing > 0

            % Create a sequence
            Seq(seq,:,blk) = GenRandSeq(L, pA(blk));

            % Check if all the patterns in this sequence are available
            if nargin == 4
                [~, nLL] = GetPatLL(nPat, Seq(seq,:,blk), ones(1,L), 0, 1, 2);
                nMissing = sum(nLL(:) == 0);
            elseif nargin < 4, nMissing = 0;
            end
        end
    end
    fprintf([bk, 'done.\n']);
end
fprintf('\n');

end
