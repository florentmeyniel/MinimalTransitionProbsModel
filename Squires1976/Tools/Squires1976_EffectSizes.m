% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% Clear the place
clear; close('all');

% For each data file
files = {'DataSquiresScience1976.mat', 'DataKolossaFIHN2013.mat'};
B = NaN(4,2);
for i = 1:2
    
    % Load data
    [Data, N, pA] = Squires1976_PrepareData(files{i});
    
    % Keep only the longest patterns
    Data = Data(end,:,:);
    Data = Data(:);
    
    % Get the patterns' list
    pats = AllSeqPattern(N);
    pats = pats{end};
    
    % Get the predictors
    P = repmat(pA, size(pats,1), 1);    % global probability
    L = 1 - mean(pats-1,2);             % local probability
    A = sum(abs(diff(pats, [], 2)), 2); % local alternations
    X = [P(:), repmat(L, numel(pA), 1), repmat(A, numel(pA), 1)];
    X = [X, ones(size(X,1), 1)];
    
    % Perform the regression
    B(:,i) = regress(Data, X);
end

% Plot the results
figure; lw = 1; fs = 20;
bar(B(1:end-1,:), 'LineWidth', lw); hold('on');
plot([0,4], zeros(1,2), 'k-', 'LineWidth', lw);
legend({'Squires et al. (1976)', 'Kolossa et al. (2013)'}, 'Location', 'SouthEast', 'Box', 'Off');
set(gca, 'XTick', 1:(size(X,2)-1), 'XTickLabel', {'Global probabilities', ...
    'Local probabilities', 'Local alternations'}, 'XTickLabelRotation', 15);
set(gca, 'Box', 'Off', 'TickDir', 'Out', 'Layer', 'Top', 'LineWidth', lw, 'FontSize', fs);
xlabel('Regressor'); ylabel('Estimated coefficient (\beta, a.u.)');
