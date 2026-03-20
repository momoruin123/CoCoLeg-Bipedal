clear; clc;
dt = 3e-2;
% load data\31_01_WS_u_squared_CoT_v_avg_BF_min_4.mat
load data\19_02_PC_u_squared_CoT_v_avg_BF_min_2.mat
configfile_BF_min;

%%
fprintf('===== Animation player =====\n');
fprintf('Intruction：\n');
% fprintf('1. select model：m / c (warmStart/Continuous)\n');
fprintf('2. Type：k_h,k_k（e.g.: 10,20）\n');
fprintf('3. Quit: q / Q \n');
fprintf('============================\n\n');

%%
while true
    % 1. Wait for Input
    user_input = input('Type stiffness（or q/Q to quit): ', 's');
    
    % 2. if quit or not
    if strcmpi(user_input, 'q')
        fprintf('The program has quitted！\n');
        break;
    end
    
    % 3. find coordinate
    try
        coord_str = strsplit(user_input, ',');
        if length(coord_str) ~= 2
            error('Coordinate format error，please type ''x,y''！');
        end
        x = str2double(coord_str{1});
        y = str2double(coord_str{2});
        current_coord = [x, y];

        idx = find(all(gridPts == current_coord, 2));    
        trajPlot = interpolate_Z2traj(config, z_final(idx,:)', [], 0, dt);
        traj_full = trajToFulltraj(config, trajPlot, 1);

        getAnimationRABBIT(traj_full, 0, config, []);

    catch ME
        fprintf('Illegal input：%s，pls retype！\n', ME.message);
        continue;
    end
end

%% Save Animation
for i = 1:15
    close all;
    trajPlot(i) = interpolate_Z2traj(config, z_final(i,:)', [], 0, dt);
    traj_full = trajToFulltraj(config, trajPlot(i), 1);
    time = datetime('now', 'Format','uuuuMMdd_HHmmss');
    time = char(time);
    gif_name = [time,'_',i];
    getAnimationRABBIT(traj_full, 0, config, []);
end
