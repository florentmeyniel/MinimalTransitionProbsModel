% This script is a toy example for Ideal Observers that estimate transition
% probabilities.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
clear; close('all');
try cd('MarkovReview'); catch, end;
addpath('Tools');
addpath('IdealObserversCode');

%% DEFINE A SEQUENCE BASED ON TRANSITION PROBABILITIES, WITH JUMPS
%  ===============================================================

% Set parameters
p     = 0.25;
p1g2  = 1-[ p  1-p  p  1-p];
p2g1  =   [ p  1-p  p  1-p];
L     =   [100 200 100 200]; % length of each chunk
N     = sum(L);
chklm = cumsum(L);
pJump = (length(p1g2) - 1) / N;

% Generate the sequence
[s, ~, gen_p1g2, gen_p2g1] = GenRandSeq(L, [p1g2', p2g1']);

%% EXAMPLE 1: COMPUTE THE IDEAL OBSERVER WITH JUMPS (HMM)
%  ======================================================

% Set parameters
in.s            = s;                % sequence
in.learned      = 'transition';     % estimate transition
in.jump         = 1;                % estimate with jumps
in.mode         = 'HMM';            % use the HMM (not sampling) algorithm
in.opt.pJ       = pJump;            % a priori probability that a jump occur at each outcome
n               = 50;               % resolution of the univariate probability grid
in.opt.pgrid    = linspace(0,1,n);  % estimation probability grid
in.opt.Alpha0   = ones(n)/(n^2);    % uniform prior on transition probabilities
in.verbose      = 1;                % to check that no default values are used.

% Compute the observer
out.HMM = IdealObserver(in);

%% EXAMPLE 2: COMPUTE THE IDEAL OBSERVER WITHOUT JUMPS
%  ===================================================

% Set parameters
clear('in');
in.s                = s;                % sequence
in.learned          = 'transition';     % estimate transition
in.jump             = 0;                % estimate with jumps
in.opt.MemParam     = {'Decay', 16};    % memory limit
in.opt.AboutFirst   = 'WithoutFirst';   % discard 1st observation for analytical solutions
n                   = 20;               % resolution of the univariate probability grid
in.opt.pgrid        = linspace(0,1,n);  % grid to return full distributions
in.opt.ReturnDist   = 1;                % Return full posterior distributions
in.opt.priorp1g2    = [1 1];            % uniform Beta prior
in.opt.priorp2g1    = [1 1];            % uniform Beta prior
in.verbose          = 1;                % to check that no default values are used.

% Compute the observer
out.FIX = IdealObserver(in);

%% PLOT THE RESULT
%  ==============

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
    subplot(7,1,1); lgd = [];
    lgd(1,:) = plot(x, gen_p1g2, '-', 'Color', col(1,:), 'LineWidth', 1.5); hold('on');
    lgd(2,:) = plot(x, gen_p2g1, '-', 'Color', col(2,:), 'LineWidth', 1.5);
    axis([1, N, 0, 1]);
    SuperposeSequence(s, 10, col(1:2,:));
    for j = 1:numel(L), plot(repmat(chklm(j),1,2), ylim, 'k:', 'LineWidth', 1); end
    legend(lgd, {'p(1|2)', 'p(2|1)'});
    ylabel({'Generative', 'probabilities and', 'simulated sequence'});

    % Full posterior
    var = {'p1g2', 'p2g1'};
    for v = 1:numel(var)
        subplot(7,1,v+1);
        d = sprintf('%s_dist', var{v});
        imagesc(x, in.opt.pgrid, out.(f{i}).(d)); hold('on');
        set(gca, 'Ydir', 'Normal');
        axis([1, N, 0, 1]);
        for j = 1:numel(L), plot(repmat(chklm(j),1,2), ylim, 'w:', 'LineWidth', 1); end
        ylabel({'Marginal distri-', sprintf('bution of \theta_{%s|%s}', var{v}(2), var{v}(4))});
    end
    
    % Mean
    subplot(7,1,4);
    var = {'p1g2', 'p2g1'};
    for v = 1:numel(var)
        m = sprintf('%s_mean', var{v});
        d = sprintf('%s_sd',   var{v});
        g = sprintf('gen_%s',  var{v});    
        SEMup = out.(f{i}).(m) + out.(f{i}).(d);
        SEMdw = out.(f{i}).(m) - out.(f{i}).(d);
        plot(x, eval(g), '--', 'Color', col(v,:), 'LineWidth', 1); hold('on');
        fill([x, flipud(x')'], [SEMup, flipud(SEMdw')'], 'k', ...
            'EdgeColor', col(v,:), 'LineWidth', 1, 'LineStyle', 'none', ...
            'FaceColor', col(v,:), 'FaceAlpha', alpha);
        plot(x, out.(f{i}).(m), '-', 'Color', col(v,:), 'LineWidth', 1.5);
    end
    axis([1, N, 0, 1]);
    for j = 1:numel(L), plot(repmat(chklm(j),1,2), ylim, 'k:', 'LineWidth', 1); end
    ylabel({'Estimation', '\theta'});

    % Log-precision
    subplot(7,1,5);
    plot(x, -log(out.(f{i}).p1g2_sd), '-', 'Color', col(1,:), 'LineWidth', 1.5); hold('on');
    plot(x, -log(out.(f{i}).p2g1_sd), '-', 'Color', col(2,:), 'LineWidth', 1.5);
    axis([1, N, 1, 4]);
    for j = 1:numel(L), plot(repmat(chklm(j),1,2), ylim, 'k:', 'LineWidth', 1); end
    ylabel({'Log-precision', '-log(\sigma_{\theta})'});

    % Surprise
    subplot(7,1,6);
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
