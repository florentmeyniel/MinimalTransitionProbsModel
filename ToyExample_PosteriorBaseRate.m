% This script is a toy example for Ideal Observers that estimate
% item frequencies.  
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
clear; close('all');
try cd('MarkovReview'); catch, end;
addpath('Tools');
addpath('IdealObserversCode');

%% DEFINE A BERNOULLI SEQUENCE WITH JUMPS
%  ======================================

% Sequence specifications
p     = 0.25;
pA    = [p, 1-p, p, 1-p, p];
L     = repmat(200, 1, numel(pA)); % length of each chunk
N     = sum(L);
chklm = cumsum(L);
pJump = (length(pA)-1) / sum(L);

% Generate the sequence
[s, gen_p1] = GenRandSeq(L, pA');

%% EXAMPLE 1: COMPUTE THE IDEAL OBSERVER WITH JUMPS (HMM)
%  ======================================================

% Set parameters
in.s            = s;                % sequence
in.learned      = 'frequency';      % estimate item frequency
in.jump         = 1;                % estimate with jumps
in.mode         = 'HMM';            % use the HMM (not sampling) algorithm
in.opt.pJ       = pJump;            % a priori probability that a jump occur at each outcome
n               = 100;              % resolution of the univariate probability grid
in.opt.pgrid    = linspace(0,1,n);  % estimation probability grid
in.opt.Alpha0   = ones(n,1)/n;      % uniform prior on transition probabilities
in.verbose      = 1;                % to check that no default values are used.

% Compute the observer
out.HMM = IdealObserver(in);

%% EXAMPLE 2: COMPUTE THE IDEAL OBSERVER WITHOUT JUMPS
%  ===================================================

% Set parameters
clear('in');
in.s                = s;                % sequence
in.learned          = 'frequency';      % estimate item frequency
in.jump             = 0;                % estimate with jumps
in.opt.MemParam     = {'Limited', 16};  % memory limit
n                   = 100;              % resolution of the univariate probability grid
in.opt.pgrid        = linspace(0,1,n);  % grid to return full distributions
in.opt.ReturnDist   = 1;                % Return full posterior distributions
in.opt.priorp1      = [1 1];            % uniform Beta prior
in.verbose          = 1;                % to check that no default values are used.

% Compute the observer
out.FIX = IdealObserver(in);

%% PLOT THE RESULT
%  ===============

% Figure options
alpha = 0.2;
x = 1:N;
col = lines(5); col = col([1,2,4,5],:);
pos = [.1 .1 .4 .75; .5 .1 .4 .75];
f = fieldnames(out);
for i = 1:numel(f)
    figure; clf; set(gcf, 'Name', out.(f{i}).FunctionName);
    set(gcf, 'Color', ones(1,3), 'Units', 'Normalized', 'Position', pos(i,:));

    % Sequence
    subplot(6,1,1); lgd = [];
    lgd(1,:) = plot(x, gen_p1, '-', 'Color', col(1,:), 'LineWidth', 1.5); hold('on');
    axis([1, N, 0, 1]);
    SuperposeSequence(s, 10, col(1:2,:));
    for j = 1:numel(L), plot(repmat(chklm(j),1,2), ylim, 'k:', 'LineWidth', 1); end
    legend(lgd, 'p(1)');
    ylabel({'Generative', 'probabilities and', 'simulated sequence'});

    % Full posterior
    subplot(6,1,2);
    imagesc(x, in.opt.pgrid, out.(f{i}).p1_dist); hold('on');
    set(gca, 'Ydir', 'Normal');
    axis([1, N, 0, 1]);
    for j = 1:numel(L), plot(repmat(chklm(j),1,2), ylim, 'w:', 'LineWidth', 1); end
    ylabel({'Marginal distri-', 'bution of \theta_{1}'});

    % Mean
    subplot(6,1,3);
    plot(x, gen_p1, '--', 'Color', col(1,:), 'LineWidth', 1); hold('on');
    SEMup = out.(f{i}).p1_mean + out.(f{i}).p1_sd;
    SEMdw = out.(f{i}).p1_mean - out.(f{i}).p1_sd;
    fill([x, flipud(x')'], [SEMup, flipud(SEMdw')'], 'k', ...
        'EdgeColor', col(1,:), 'LineWidth', 1, 'LineStyle', 'none', ...
        'FaceColor', col(1,:), 'FaceAlpha', alpha);
    plot(x, out.(f{i}).p1_mean, '-', 'Color', col(1,:), 'LineWidth', 1.5); hold('on');
    axis([1, N, 0, 1]);
    for j = 1:numel(L), plot(repmat(chklm(j),1,2), ylim, 'k:', 'LineWidth', 1); end
    ylabel({'Estimation', '\theta'});

    % Log-precision
    subplot(6,1,4);
    plot(x, -log(out.(f{i}).p1_sd), '-', 'Color', col(1,:), 'LineWidth', 1.5); hold('on');
    axis([1, N, 1, 4]);
    for j = 1:numel(L), plot(repmat(chklm(j),1,2), ylim, 'k:', 'LineWidth', 1); end
    ylabel({'Log-precision', '-log(\sigma_{\theta})'});

    % Surprise
    subplot(6,1,5);
    plot(x, out.(f{i}).surprise, '.', 'Color', col(3,:), 'MarkerSize', 10); hold('on');
    axis([1, N, 0, 4]);
    for j = 1:numel(L), plot(repmat(chklm(j),1,2), ylim, 'k:', 'LineWidth', 1); end
    ylabel({'Surprise', '-log(\theta)'});

    % Update
    subplot(7,1,7);
    plot(x, out.(f{i}).distUpdate, '-', 'Color', col(4,:), 'LineWidth', 1.5); hold('on');
    axis([1, N, 0, 1]);
    for j = 1:numel(L), plot(repmat(chklm(j),1,2), ylim, 'k:', 'LineWidth', 1); end
    xlabel('Observation #'); ylabel({'Update', 'D_{KL}(\theta_{k-1}||\theta_{k})'});
end
