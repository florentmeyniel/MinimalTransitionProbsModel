function p = SuperposeSequence( s, ptsize, col )
%SUPERPOSESEQUENCE allows to overlap a binary sequence over an existing
%plot.
%   - "s": a 1xN vector specifying the binary sequence to plot.
%   - "ptsize": a scalar sprcifying the size of the dots.
%   - "col": a 2x3 matrix specifying the RGB color vector for each stimulus
%   category.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% Fill the empty inputs
if nargin < 2, ptsize = 10;    end
if nargin < 3, col = lines(2); end

% Get current figure
lim = get(gca, 'ylim');
figure(gcf); hold('on');

% Plot the 
p = [];
for i = 1:2
    try p(i,:) = plot(find(s == i), lim(i), '.', 'Color', col(i,:), 'MarkerSize', ptsize); catch, end;
end

% Scale the plot
ylim(lim);
set(gca, 'Layer', 'Bottom');

end
