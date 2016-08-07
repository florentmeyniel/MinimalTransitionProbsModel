function lgd = PlotTree( m, N, col, sd, n, stim, fs, lw, xoffset, letfun, gap )
%PLOTTREE plots a matrix as a tree (with overlapping labels).
%   - "m": a NxS matrix (where S = size(ff2n(N),1)/2) specifying the
%   averaged likelihood value associated to each pattern. These are the
%   values that will shape the tree. This is the only mandatory input
%   variable.
%   - "N": the maximum length of the patterns to plot.
%   - "col": the color of the tree to plot.
%   - "sd": a NxS matrix specifying the likelihood standard deviations.
%   - "n": a NxS matrix specifying the number of observations used to
%   generate "m" (the last 2 inputs are used to plod error bars on the
%   figure).
%   - "stim": the letter used to make the labels (1: A, 2: B, 3: X).
%   - "fs": the size of the labels' font.
%   - "lw": the tickness of the lines.
%   - "xoffset": the x-tick at which the tree must start (the tree will be
%   plotted at the left of this point).
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

%% INITIALIZATION
%  ==============

% Default options
if nargin < 11, gap = 0.1;
    if nargin < 10, letfun = @upper;
        if nargin < 9, xoffset = 0;
            if nargin < 8, lw = 2;
                if nargin < 7, fs = 12;
                    if nargin < 6, stim = 1;
                        if nargin < 5, n = [];
                            if nargin < 4, sd = [];
                                if nargin < 3, col = lines(1);
                                    if nargin < 2, N = size(m,1); end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% Abscissa vector
x = (N-1+xoffset):-1:xoffset;

% If the empty cells have been filled with zeros, replace them by NaNs
if ~isnan(sum(m(:)))
    for i = 1:size(m,1)
        ii = size(ff2n(i),1)/2;
        m(i,ii+1:end) = NaN;
    end
end

% Define the orientation of the text
if     gap <  0, orient = 'right';
elseif gap >= 0, orient = 'left'; end

%% PLOT THE TREE
%  =============

% For each
for k = 1:N
    if k < N
        % Get start (k) point (shorter sequence) and end (k+1) point
        % (longer sequence)
        startpoint  = m(k  , ~isnan(m(k  ,:)));
        destination = m(k+1, ~isnan(m(k+1,:)));

        % The 1st part (up to end-1) of the longer sequence is identical to
        % the second. This is because the longer sequence is derived from
        % the shorter by adding as a 1st element A (1st half) or B (second
        % half).
        indhalf = length(destination) / 2;
    end
    for iStart = 1:length(startpoint)
        
        % Plot the dots
        if k < N
            % Longer sequences starting with A
            lgd = plot([x(k), x(k+1)], [m(k,iStart), m(k+1,iStart)], ...
                '.-', 'Color', col, 'LineWidth', lw, ...
                'MarkerFaceColor', col, 'MarkerSize', lw*15); hold('on');

            % Longer sequences starting with B
            plot([x(k), x(k+1)], [m(k,iStart), m(k+1,indhalf+iStart)], ...
                '.-', 'Color', col, 'LineWidth', lw, ...
                'MarkerFaceColor', col, 'MarkerSize', lw*15); hold('on');
        end
        
        % Plot the error bars
        if ~isempty(sd) && ~isempty(n)
            errorbar(x(k), m(k,iStart), sd(k,iStart)./sqrt(n(k,iStart)), ...
                'Color', col, 'LineWidth', lw); hold('on');
        end
    end
end

%% OVERLAP PATTERNS' LABELS
%  ========================

% Text labels
if     stim == 1, stim = {'A', 'B'};
elseif stim == 2, stim = {'B', 'A'};
elseif stim == 3, stim = {'X', 'Y'};
elseif stim == 4, stim = {'X', 'O'};
end
stim = letfun(stim);

% Get all combinations
seq = AllSeqPattern(N);

% Alone stimulus
text(x(1)+gap, m(1,1), sprintf(strcat(stim{seq{1}})),...
    'HorizontalAlignment', orient, 'FontSize', fs); hold('on');

% Combined stimuli
for k = 1:N
    if k < N
        startpoint  = m(k  , ~isnan(m(k  ,:)));
        destination = m(k+1, ~isnan(m(k+1,:)));
        indhalf = length(destination) / 2;
    end
    for iStart = 1:length(startpoint)
        if k < N
            text(x(k+1)+gap, m(k+1,iStart), ...
                sprintf(strcat(stim{seq{k+1}(iStart,:)})), ...
                'HorizontalAlignment', orient, 'FontSize', fs); hold('on');
            text(x(k+1)+gap, m(k+1,indhalf+iStart), ...
                sprintf(strcat(stim{seq{k+1}(indhalf+iStart,:)})), ...
                'HorizontalAlignment', orient, 'FontSize', fs); hold('on');
        end
    end
end

%% PLOT OPTIONS
%  ============

xlim([x(end) - 1, x(1) + 1]);
set(gca, 'XTick', fliplr(x));%, 'XTickLabel', num2str((N:-1:1)'));

end
