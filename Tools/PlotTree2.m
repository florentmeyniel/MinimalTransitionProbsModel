function lgd = PlotTree2( m, N, col, lwd, xoffset, pt, dec, treecol, lwt )
%PLOTTREE plots a matrix as a tree (with small flages).
%   - "m": a NxS matrix (where S = size(ff2n(N),1)/2) specifying the
%   averaged likelihood value associated to each pattern. These are the
%   values that will shape the tree. This is the only mandatory input
%   variable.
%   - "N": the maximum length of the patterns to plot.
%   - "col": the color of the ____.
%   - "fs": the size of the labels' font.
%   - "lw": the tickness of the lines.
%   - "xoffset": the x-tick at which the tree must start (the tree will be
%   plotted at the left of this point).
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

if nargin < 9, lwt = 1;
    if nargin < 8, treecol = repmat(128,1,3)./255;
        if nargin < 7, dec = 0.1;
            if nargin < 6, pt = 7;
                if nargin < 5, xoffset = 0;
                    if nargin < 4, lwd = 0.5;
                        if nargin < 3, col = [255 000 000; 000 000 255] ./ 255;
                            if nargin < 2, N = size(m,1);
                            end
                        end
                    end
                end
            end
        end
    end
end

% X coordinates
xmed = (N-1+xoffset):-0.5:xoffset;

% Get all combinations
seq = AllSeqPattern(N);

% For each pattern length
for i = 1:N
    
    % Draw lines
    if i < N
        startpoint  = m(i  , ~isnan(m(i  ,:)));
        destination = m(i+1, ~isnan(m(i+1,:)));
        indhalf = length(destination) / 2;
        
        for iStart = 1:length(startpoint)
            lgd = plot([xmed(i), xmed(i+1)], [m(i,iStart), m(i+1,iStart)],   '-', 'Color', treecol, 'LineWidth', lwt); hold('on');
            plot([xmed(i), xmed(i+1)], [m(i,iStart), m(i+1,indhalf+iStart)], '-', 'Color', treecol, 'LineWidth', lwt);
        end
    end
    
    % For each pattern
    for j = 1:size(seq{i},1)
        pattern = seq{i}(j,:);
        n = numel(pattern);
        
        % Define positions of the circles centered on the common x position
        if i > 1
            xarray = (dec:dec:(n*dec));
            xarray = xarray - mean(xarray);
        else xarray = 0;
        end
        
        % For each stimulus
        for k = 1:n
            
            % Display a colored circle
            plot(xmed(i)+xarray(k), m(i,j), 'o', 'MarkerSize', pt, 'LineWidth', lwd, ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', col(pattern(k),:)); hold('on');
        end
    end
end

% Scale the tree
xlim([xmed(end) - 1, xmed(1) + 1]);
set(gca, 'XTick', fliplr(xmed));

end
