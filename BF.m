clear;clc;

% configfile_BF_m;
% configfile_BF_min;
configfile_BF_min_flight;
% configfile_BF_min_2;
dt = 3e-3;
traj = flow(config, dt);
traj_full = trajToFulltraj(config,traj);
getAnimationRABBIT(traj_full, 0, config,[]);
