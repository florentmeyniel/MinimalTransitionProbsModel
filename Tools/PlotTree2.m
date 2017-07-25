function lgd = PlotTree2( m, N, col, lwd, xoffset, pt, dec, treecol, lwt, ses )
%PLOTTREE plots a matrix as a tree (with small flags).
%   - "m": a NxS matrix (where S = size(ff2n(N),1)/2) specifying the
%   (averaged) likelihood value or surprise level associated to each
%   pattern. These are the values that will shape the tree. N.B. This is
%   the only mandatory input variable.
%   - "N": the maximum length of the patterns to plot.
%   - "col": a 2x3 color matrix specifying the colors 
%   - "lwd": width of the line surrounding the circles composing the flags.
%   - "xoffset": the x-tick at which the root of the tree should be set
%   N.B. the tree will be plotted at the left of this point.
%   - "pt": size of the circles composing the flags.
%   - "dec": spacing between 2 circles in the flags.
%   - "treecol": the color of the tree's branches.
%   - "lwt": the tickness of the lines.
%   - "ses": whether or not (boolean) underlie the last stimulus of each
%   pattern (such that it is made explicit that the y position of each flag
%   represents the measure (specified in m) corresponding to the last
%   element of each pattern.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% Complete inputs
if nargin < 10, ses = false;
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
end
warning('OFF');

% X-axis coordinates
xmed = (N-1+xoffset):-0.5:xoffset;

% Get all combinations
seq = AllSeqPattern(N);

% For each pattern length
for i = 1:N
    
    % Draw lines
    % ~~~~~~~~~~
    if i < N
        startpoint  = m(i  , ~isnan(m(i  ,:)));
        destination = m(i+1, ~isnan(m(i+1,:)));
        indhalf = length(destination) / 2;
        for iStart = 1:length(startpoint)
            lgd = plot([xmed(i), xmed(i+1)], [m(i,iStart), m(i+1,iStart)], ...
                '-', 'Color', treecol, 'LineWidth', lwt); hold('on');
            plot([xmed(i), xmed(i+1)], [m(i,iStart), m(i+1,indhalf+iStart)], ...
                '-', 'Color', treecol, 'LineWidth', lwt);
        end
    end
    
    % Draw flags
    % ~~~~~~~~~~
    for j = 1:size(seq{i},1)
        pattern = seq{i}(j,:);
        n = numel(pattern);
        
        % Each stimulus is displayed as a colored circle
        if i > 1
            xarray = (dec:dec:(n*dec)); % x positions of each circle
            xarray = xarray - mean(xarray);
        else, xarray = 0;
        end
        for k = 1:n
            plot(xmed(i)+xarray(k), m(i,j), 'o', 'MarkerSize', pt, 'LineWidth', lwd, ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', col(pattern(k),:)); hold('on');
        end
        
        % Underlie the last stimulus
        if ses
            text(xmed(i)+xarray(n), m(i,j), '_', ...
                'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        end
    end
end

% Scale the tree
xlim([xmed(end)-1, xmed(1)+1]);
set(gca, 'XTick', fliplr(xmed));

end
