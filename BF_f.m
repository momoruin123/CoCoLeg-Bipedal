clear all;clc;
configfile_BF_floating;
dt = 3e-2;
traj = flow(config, dt);

getAnimationRABBIT(traj, 0, config,[]);