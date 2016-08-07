function Y = Squires1976_KolossaLLpat(s)
% return the surprise at observing "a" as the last stimulus.
% ex: LLpat('bbba')
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% check that the last, observed stimulus is a
if ~strcmp(s(end)', 'a')
    error('the last stimulus should be a')
end

% convert the sequence into 0 and 1
s_num = zeros(1,length(s));
s_num(strfind(s, 'a')) = 1;
s = s_num;

% remove the last observation, since it should not be taken into account in
% the inference of the hidden statistics.
s = s(1:end-1);

if length(s) == 3
    % complete so the there are 4 previous observations
    p_seq1 = Squires1976_KolossaModel([0 s]);
    p_seq2 = Squires1976_KolossaModel([1 s]);
    Y = -0.5*(log2(p_seq1(end))+log2(p_seq2(end)));
    
elseif length(s) == 2
    % complete so the there are 4 previous observations
    p_seq1 = Squires1976_KolossaModel([0 0 s]);
    p_seq2 = Squires1976_KolossaModel([1 0 s]);
    p_seq3 = Squires1976_KolossaModel([0 1 s]);
    p_seq4 = Squires1976_KolossaModel([1 1 s]);
    Y = -0.25*(log2(p_seq1(end))+log2(p_seq2(end))+log2(p_seq3(end))+log2(p_seq4(end)));
elseif length(s) == 1
    % complete so the there are 4 previous observations
    p_seq1 = Squires1976_KolossaModel([0 0 0 s]);
    p_seq2 = Squires1976_KolossaModel([0 1 0 s]);
    p_seq3 = Squires1976_KolossaModel([0 0 1 s]);
    p_seq4 = Squires1976_KolossaModel([0 1 1 s]);
    
    p_seq5 = Squires1976_KolossaModel([1 0 0 s]);
    p_seq6 = Squires1976_KolossaModel([1 1 0 s]);
    p_seq7 = Squires1976_KolossaModel([1 0 1 s]);
    p_seq8 = Squires1976_KolossaModel([1 1 1 s]);
    
    Y = -0.125*(log2(p_seq1(end))+log2(p_seq2(end))+log2(p_seq3(end))+log2(p_seq4(end))+...
                log2(p_seq5(end))+log2(p_seq6(end))+log2(p_seq7(end))+log2(p_seq8(end)));
    
end

end
