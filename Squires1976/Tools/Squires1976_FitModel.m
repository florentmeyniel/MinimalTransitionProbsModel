function [ MSE_nW, MSE_W, R2_nW, R2_W, sc_nW, sc_W, beta_nW, beta_W ] = ...
    Squires1976_FitModel( X, Y, W )
%SQUIRES1976_FITMODEL performs weighted and not weighted linear regression.
%   - "X": the predictor matrix.
%   - "Y": the data matrix.
%   - "W": the weight matrix.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 

% Complete the inputs
I = ~isnan(X);
if nargin < 3
    W = NaN(size(X));
    W(I) = 1;
end

% Check that the inputs have the same size
sz = [size(X); size(Y); size(W)];
if sum(diff(sz,2)) ~= 0
    error('The inputs do not have the same size');
end

% Deal with NaN values
X = X(I);
Y = Y(I);
W = W(I);

% Deal with infine values
X = X(~isinf(X));
Y = Y(~isinf(X));
W = W(~isinf(X));

% Get weight on data across tested probabilities
W = (W * ones(1,numel(X))) .* eye(numel(X)) .^ (1/1);

% Add an offset to X
X = [X, ones(numel(X),1)];

% Estimate coefficients
beta_nW = (X'     * X) \ X'     * Y;                % not weighted
beta_W  = (X' * W * X) \ X' * W * Y;                % weighted

% Generate predictions
yhat_nW = X * beta_nW;                              % not weighted
yhat_W  = X * beta_W;                               % weighted

% Measure the error
MSE_nW = mean(( Y - yhat_nW) .^ 2);                 % not weighted
MSE_W  = mean(((Y - yhat_W ) .^ 2) .* sum(W, 2));   % weighted

% Get the proportion of explained variance
R2_nW = var(yhat_nW) / var(Y);                      % not weighted
R2_W  = var(yhat_W)  / var(Y);                      % weighted

% Export the scaled predictor in the same format as the inputs
sc_nW = NaN(size(I));
sc_nW(I) = yhat_nW;
sc_W = NaN(size(I));
sc_W(I) = yhat_W;

end
