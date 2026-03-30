clear;
clc;
%% BF_min: Hopping in place with minimum coordiantes
% configfile_BF_min;
configfile_BF_run;
config.optConfig.print_level = 5;

%% Configuration
dt = 3e-3;
config.paramValues.k_h = 20;
config.paramValues.k_k = 20;

%% Initial guess
% trajInit = flow(config, dt);

trajInit = initial_guess(config, 1, 0.5, 1);
traj_full = trajToFulltraj(config, trajInit, 1);
getAnimationRABBIT(traj_full, 0, config,[]);

% load result\trajPlot_20260211_113035.mat

%% Optimization
% Set up and solve optimization
[Z, P, w, r, constraintNames, h, bounds, boundNames] = getConstraintsAndCosts(config, trajInit);

%%
[solutionIPOPT,  stats] = solveCasadi(config, Z, P, r, h, w, bounds, trajInit);

%%
solution = full(solutionIPOPT.x);
trajPlot = interpolate_Z2traj(config, solution, [], 0, dt);
traj_full = trajToFulltraj(config, trajPlot, 1);

%%
% plot_traj(config, trajPlot, [], [], []);

%% Save Animation
save_flag = false;
file_name = 'run_gait';
% file_name = mfilename;
time = datetime('now', 'Format','uuuuMMdd_HHmmss');
time = char(time);
gif_name = [file_name,'_',time];
if save_flag
    getAnimationRABBIT(traj_full, 0, config, gif_name);
    save(['trajPlot','_',time], 'trajPlot');
    save(['solution','_',time], 'solution');
    save(['traj_full','_',time], 'traj_full');
else
    getAnimationRABBIT(traj_full, 0, config,[]);
end
