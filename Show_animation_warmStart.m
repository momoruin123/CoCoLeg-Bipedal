clear; clc;
dt = 3e-2;
% load data\31_01_WS_u_squared_CoT_v_avg_BF_min_4.mat
% load data\12_03_WS_u_squared_CoT_v_avg_BF_run_2.mat
load data\31_03_WS_u_squared_CoT_v_avg_BF_run_1.mat

configfile_BF_run;

%% Save Animation
for i = 1:10
    close all;
    trajPlot = interpolate_Z2traj(config, Z_array(i,:)', [], 0, dt);
    traj_full = trajToFulltraj(config, trajPlot, 0);
    time = datetime('now', 'Format','uuuuMMdd_HHmmss');
    time = char(time);
    gif_name = [time,'_',i];
    getAnimationRABBIT(traj_full, 0, config, []);
end
