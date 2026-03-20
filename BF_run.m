clear;clc;

% configfile_BF_m;
configfile_BF_run;
% configfile_BF_min;
dt = 3e-3;
% traj = flow(config, dt);
% traj_full = trajToFulltraj(config, traj, 1);

trajInit = initial_guess(config, 1);
traj_full = trajToFulltraj(config, trajInit, 1);
getAnimationRABBIT(traj_full, 0, config,[]);
