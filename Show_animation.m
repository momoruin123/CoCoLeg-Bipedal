clear; clc;
dt = 3e-2;
% load data\31_01_WS_u_squared_CoT_v_avg_BF_min_4.mat
load data\11_02_WS_u_squared_CoT_v_avg_BF_min_5.mat
configfile_BF_min;

%% Save Animation
for i = 1:15
    close all;
    trajPlot(i) = interpolate_Z2traj(config, Z_array(i,:)', [], 0, dt);
    traj_full = trajToFulltraj(config, trajPlot(i), 1);
    time = datetime('now', 'Format','uuuuMMdd_HHmmss');
    time = char(time);
    gif_name = [time,'_',i];
    getAnimationRABBIT(traj_full, 0, config, []);
end
