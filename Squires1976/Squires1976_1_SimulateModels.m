% Run the simulation corresponding to Squires et al 1976 experiment.
% 
% Copyright 2016 Florent Meyniel & Maxime Maheu 
 
%% INITIALIZATION
%  ==============

% Clear the place
clear; close('all'); clc;

% Get the functions
addpath(genpath('Squires1976'));
try cd('Squires1976'); catch, end

% Load data
[~, nPat, p1] = Squires1976_PrepareData('DataSquiresScience1976.mat');
nP = numel(p1);

% Reorder the blocks
p1 = sort(p1, 2, 'descend');

%% GENERATE SEQUENCES
%  ==================

% Display what we are doing
fprintf('GENERATING SEQUENCES...\n');

% Define sequences' parameters
L    = 200;
nSeq = 200;

% Generating the stimuli sequences
Spec.StimSeq = Squires1976_GenerateSequences(p1, L, nSeq, nPat);

%% SIMULATE THE PRINCIPLED OBSERVERS
%  =================================

% Display what we are doing
fprintf('\nSIMULATING THE PRINCIPLED OBSERVERS...\n\n');

% Define the models to fit
Info.Models = {'Stimulus',    'Windowed'; ...
               'Stimulus',    'Leaky'; ...
               'Stimulus',    'Dynamic'; ...
               'Alternation', 'Windowed'; ...
               'Alternation', 'Leaky'; ...
               'Alternation', 'Dynamic'; ...
               'Transitions', 'Windowed'; ...
               'Transitions', 'Leaky'; ...
               'Transitions', 'Dynamic'; ...
               'Descriptive', 'Squires et al. (1976)'; ...
               'Descriptive', 'Kolossa et al. (2013)'};
nMod = size(Info.Models, 1);

% Define the quantity to estimate
Info.Quty  = {'Surprise', 'Update'};
nQuty = numel(Info.Quty);

% Prepare variables for the simulations' results
Spec.ParamGrid   = cell(1,nMod);
Simu = cell(1,nMod);

% Shared options among observers
n = 20;
shr                  = [];
shr.opt.pgrid        = linspace(0,1,n);
shr.verbose          = 0;
shr.opt.ReturnDist	 = 0;
shr.opt.priorp1      = ones(1,2);
shr.opt.priorp1g2    = ones(1,2);
shr.opt.priorp2g1    = ones(1,2);

% For each memory type and each type of statistic to learn
for m = 1:nMod-2
    
    % Get the model
    stat = Info.Models{m,1};
    intg = Info.Models{m,2};
    fprintf('Learning "%s" with a "%s" integration:\n', stat, intg);

    % Get the parameters' grid
    if     strcmpi(intg, 'Windowed') || strcmpi(intg, 'Leaky')
        Spec.ParamGrid{m} = [1:1:25 27 29 32 35 38 41 44 48 55 62 71 82 98 116 145 170 200];
    elseif strcmpi(intg, 'Dynamic')
        Spec.ParamGrid{m} = [0.0005 0.002 0.005 0.0083 0.013 0.019 ...
              0.025 0.033 0.042 0.052 0.063 0.074 0.087 0.101 0.116 0.132 ...
              0.149 0.167 0.186 0.207 0.228 0.250 0.273 0.300 0.323 0.349 ...
              0.377 0.405 0.434 0.465 0.496 0.529 0.563 0.597 0.633 0.669 ...
              0.707 0.746 0.786 0.826 0.911 1];
    end
    nParam = numel(Spec.ParamGrid{m});

    % Variables used for printing purpose
    np = numel(sprintf('%i', nParam));
    nd = numel(sprintf('%i', nSeq));
    rp = 11+(2*nd);
    bk = repmat('\b', [1,rp]);
    
    % Prepare the output
    Simu{m} = NaN(L, nP, nParam, nQuty, nSeq);

    % For each parameter value
    for p = 1:nParam
        
        % Get the parmater value
        param = Spec.ParamGrid{m}(p);
        fprintf([' - Parameter %', num2str(np), '.0f/%', num2str(np), ...
            '.0f (%g):\t%s'], p, nParam, Spec.ParamGrid{m}(p), repmat(' ', [1,rp]));

        % For each stimuli sequence
        for seq = 1:nSeq
            fprintf([bk, 'sequence %', num2str(nd), '.0f/%' , num2str(nd), '.0f\n'], seq, nSeq);

            % For each global probability of A
            for blk = 1:nP

                % Get the default (shared) inputs
                in = shr;

                % Get the sequence
                in.s = Spec.StimSeq(seq,:,blk);
                if strcmpi(stat, 'Alternation')
                    in.s = abs(diff(in.s));
                    in.s = abs(in.s-2);
                end
                % Note that in the case of the observer learning
                % repetitions versus alternations, the sequence is
                % transformed in the way that repetitions are coded
                % with 2s and alternations with 1s.

                % Get the statistic to learn and the prior on the
                % distribution
                if     strcmpi(stat, 'Stimulus')
                    in.learned      = 'frequency';
                    in.opt.Alpha0	= ones(n,1)/n;
                elseif strcmpi(stat, 'Alternation')
                    in.learned      = 'frequency';
                    in.opt.Alpha0	= ones(n,1)/n;
                elseif strcmpi(stat, 'Transitions')
                    in.learned      = 'transition';
                    in.opt.Alpha0	= ones(n)/(n^2);
                end

                % Get the memory type of the observer
                if     strcmpi(intg, 'Windowed')
                    in.jump         = 0;
                    in.opt.MemParam = {'Limited', param};
                elseif strcmpi(intg, 'Leaky')
                    in.jump         = 0;
                    in.opt.MemParam = {'Decay', param};
                elseif strcmpi(intg, 'Dynamic')
                    in.jump         = 1;
                    in.mode         = 'HMM';
                    in.opt.pJ       = param;
                end

                % Run the observer
                out = IdealObserver(in);

                % Get the variables we want to fit
                if strcmpi(stat, 'Alternation')
                    surprise = [NaN, out.surprise];
                    update   = [NaN, out.distUpdate];
                else
                    surprise = out.surprise;
                    update   = out.distUpdate;
                end

                % Sauver ca a la place des patterns, retrouver les
                % patterns dans les autres scripts
                Simu{m}(:,blk,p,1,seq) = surprise;
                Simu{m}(:,blk,p,2,seq) = update;
            end
        end
        fprintf([bk, 'done.\n']);
    end
    fprintf('\n');
end

%% SIMULATE THE DESCRIPTIVE OBSERVERS
%  ==================================

% Display what we are doing
fprintf('SIMULATING THE DESCRIPTIVE OBSERVERS...\n\n');

% Variables used for printing purpose
nd = numel(sprintf('%i', nSeq));
rp = 11+(2*nd);
bk = repmat('\b', [1,rp]);
    
% For each global probability of A
for blk = 1:nP
    fprintf(' - Block %i/%i [p(1) = %1.2f]:\t%s', blk, nP, p1(blk), repmat(' ', [1,rp]));
    
    % For each stimuli sequence
    for seq = 1:nSeq
        fprintf([bk, 'sequence %', num2str(nd), '.0f/%' , num2str(nd), '.0f\n'], seq, nSeq);
        
        % Model from Squires et al. (1976)
        [~, Simu{end-1}(:,blk,1,1,seq)] = Squires1976_SquiresModel(...
            Spec.StimSeq(seq,:,blk), ... % the stimuli sequence
            p1(blk)); ... % the global probability of events is magically provided

        % Model from Kolossa et al. (2013)
        [~, Simu{end}(:,blk,1,1,seq)] = Squires1976_KolossaModel(...
            Spec.StimSeq(seq,:,blk)); % the stimuli sequence
    end
    fprintf([bk, 'done.\n']);
end

%% SAVE THE RESULTS
%  ================

save('Squires1976_SimulateModels.mat', '-v7.3', 'Info', 'Spec', 'Simu');
