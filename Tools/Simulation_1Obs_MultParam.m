function Simulation_1Obs_MultParam(jobfile)
% This function execute the estimation by one observer type, and loop over
% different values of the leak parameters. The parameters of the observer
% and the grid of leak values is spefied in a mat file that is loaded (and
% whose name serve as an input to this function).
%
% Copyright 2016 Florent Meyniel & Maxime Maheu

fprintf('running: %s', jobfile)

% add path 
addpath ../IdealObserversCode/
addpath ../Tools/

% load job, including parameters of the observer
load(jobfile)

% run the ideal observer estimation
n_param = numel(grid_leak);
ObsLL = nan(n_param, numel(in.s));
for i_param = 1:n_param
    fprintf('\n %d / %d', i_param, n_param)
    in.opt.MemParam     = {'Decay', grid_leak(i_param)};  % use a limited memory (exponential forgetting)
    out                 = IdealObserver(in);
    ObsLL(i_param,:)    = out.surprise;
end

% save the results
save([jobfile(1:end-10), '_res.mat'], 'ObsLL')