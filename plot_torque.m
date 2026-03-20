clear; clc;
import BF_min.SingleStance.*

%% Load data set
load data\31_01_WS_u_squared_CoT_v_avg_BF_min_4.mat
configfile_BF_min;

%% Model parameter configuration
[~, nq_m, nu, ~] = get_sizes();
dt = 3e-2;
params = getModelParameters(config, 0, 0);
paramName = {'m', 'g'};
p_idx = getParameterIndex(paramName, config);
m = params(p_idx(1));
g = params(p_idx(2));

%% Tau_spring cauculation and plotting
% initialization
trajPlot = cell(15,1);
tau_spring = cell(15,1);
curve_names = {'alpha_l','alpha_r','beta_l','beta_r','alpha_l', 'beta_l'};
v_avg = gridPoint_array(:,1);
% save path
current_path = fileparts(mfilename('fullpath'));
save_folder = fullfile(current_path, 'temp');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

for i = 1:15
    trajPlot{i} = interpolate_Z2traj(config, Z_array(i,:)', [], 0, dt);

    % Extract 
    traj_vector = trajToVector(config, trajPlot{i});
    x = traj_vector.x;
    u = traj_vector.u;
    t = traj_vector.t;
    q_act = x(:,2:5);
    dq_act = x(:, 7:10);
    p = gridPoint_array(:,2:3)*m*g;
    params(end-1:end) = p(i,:)';
    num_x = size(x,1);
    tau_spring_temp = zeros(num_x, 5);

    % Spring forces
    for j = 1:num_x
        x_j = x(j,:)';
        tau_spring_temp(j,:) = tau_Spring_m(x_j, params);
    end
    tau_spring{i} = tau_spring_temp;

    % plot figure
    fig = figure('Position',[100, 100, 1200, 800]);
    set(fig,'Name',sprintf('v_avg = %.2f',v_avg(i)),'NumberTitle','off');

    for j = 1:6
        ax = subplot(3,2,j);
        hold on; grid on; box on;
        set(ax,'FontSize',10);

        if j <=4
            % Extract data
            u_j = u(:,j);
            tau_spring_j = tau_spring{i}(:,j+1);
            t_j = t;

            % plot curves
            yyaxis left;
            plot(t_j, tau_spring_j, 'LineWidth',1.5, 'DisplayName',sprintf('tau_{spring}'));
            ylabel('tau_{spring} (N·m)','FontSize',10);
            y1 = ylim; % 获取左侧轴原始范围 [y1_min, y1_max]


            yyaxis right;
            plot(t_j, u_j, 'LineWidth',1.5, 'DisplayName',sprintf('u'));
            ylabel('u (N・m)','FontSize',10);
            y2 = ylim; % 获取右侧轴原始范围 [y2_min, y2_max]

            % 计算以0为基准的对称比例，强制0刻度对齐
            y1_abs = max(abs(y1)); % 左侧轴离0的最大绝对值
            y2_abs = max(abs(y2)); % 右侧轴离0的最大绝对值

            % 【修正语法】先激活轴，再设ylim，所有MATLAB版本兼容
            yyaxis left;
            ylim([-y1_abs, y1_abs]);  % 左侧轴设为对称范围，0在正中间
            yyaxis right;
            ylim([-y2_abs, y2_abs]); % 右侧轴设为对称范围，0在正中间
        else
            idx = (j-5)*2+1;
            q_act_j = q_act(:,idx);
            plot(t_j(1:51), q_act_j(1:51), 'LineWidth',1.5, 'color', [1,0,0], 'DisplayName',sprintf('%s', curve_names{j}));
            plot(t_j(52:end), q_act_j(52:end), 'LineWidth',1.5, 'color', [0,0,1], 'DisplayName',sprintf('%s', curve_names{j}));
            ylabel(sprintf('%s (rad)', curve_names{j}),'FontSize',10);
            y1 = ylim;
            y1_abs = max(abs(y1));
            ylim([-y1_abs, y1_abs]);
        end
        % subtitles
        title(curve_names{j},'FontSize',10,'FontWeight','bold');
        xlabel('time (s)','FontSize',10);
        if j == 1
            legend('Location','northwest','FontSize',9); 
        end
    end
    % title
    sgtitle(sprintf('v_{avg} = %.2f m/s k_h=%.2f k_k=%.2f',v_avg(i),p(i,1),p(i,2)),'FontSize',18,'FontWeight','bold','Color','k');
  
    % save figure
    save_name = sprintf('TorqueComparison_vavg%.2f.png', v_avg(i));
    save_path = fullfile(save_folder, save_name);
    exportgraphics(fig, save_path, 'Resolution', 300);

    close;
end

