%% Loading 
configfile_BF_min;

load trajPlot_20260120_232037.mat

%% config
nq = 5;
nx = 10;
nu = 4;

v = 0.7;
N = 150;

%% States
close all;

x1 = trajPlot.x;
t = trajPlot.t;
name = config.defaultStateNames;
unit = {'\theta [rad]', '\omega [rad/s]'};

% Data Phase 1 
time_1 = t;

% Data Phase 2
end_t = time_1(end);
time_2 = end_t + time_1;
x2 = x1;
x2(:,2) = x1(:,3);
x2(:,3) = x1(:,2);
x2(:,4) = x1(:,5);
x2(:,5) = x1(:,4);
x2(:,7) = x1(:,8);
x2(:,8) = x1(:,7);
x2(:,9) = x1(:,10);
x2(:,10) = x1(:,9);

for i = 1:nx
    % Figure config
    figure(i);
    hold on;grid on;
    set(gcf, 'Position', [100, 100, 800, 300]);

    data1 = x1(:,i);
    data2 = x2(:,i);
    plot(time_1, data1, 'LineWidth', 2, 'DisplayName', 'On left leg')
    plot(time_2, data2, 'LineWidth', 2, 'DisplayName', 'On right leg')

    xlabel('Time [s]');
    if i > nq
        j = 2;
    else
        j = 1;
    end
    ylabel(unit{j});

    title(name{i});
    legend('show'); % 显示图例以区分颜色
    
    savefig(figure(i), sprintf('%s_v%g_N%g.fig', name{i},v,N));% save .fig

    filename_png = sprintf('%s_v%g_N%g.png', name{i},v,N);% save .png
    exportgraphics(gcf, filename_png, ...
    'Resolution', 300, ...       % 300DPI高清
    'BackgroundColor', 'white'); % 白色背景（适配论文/报告）
end

%% Input
close all;

u1 = trajPlot.u;
t = trajPlot.t;
name = config.defaultInputNames;

% Data Phase 1 
time_1 = t;

% Data Phase 2
end_t = time_1(end);
time_2 = end_t + time_1;

u2 = u1;
u2(:,1) = u1(:,2);
u2(:,2) = u1(:,1);
u2(:,3) = u1(:,4);
u2(:,4) = u1(:,3);

for i = 1:nu
    % Figure config
    figure(i);
    hold on;grid on;
    set(gcf, 'Position', [100, 100, 800, 300]);

    data1 = u1(:,i);
    data2 = u2(:,i);
    plot(time_1, data1, 'LineWidth', 2, 'DisplayName', 'On left leg')
    plot(time_2, data2, 'LineWidth', 2, 'DisplayName', 'On right leg')

    xlabel('Time [s]');
    ylabel('Torque [N·m]');

    title(name{i});
    legend('show'); % 显示图例以区分颜色
    
    savefig(figure(i), sprintf('%s_v%g_N%g.fig', name{i},v,N));% save .fig

    filename_png = sprintf('%s_v%g_N%g.png', name{i},v,N);% save .png
    exportgraphics(gcf, filename_png, ...
    'Resolution', 300, ...       % 300DPI高清
    'BackgroundColor', 'white'); % 白色背景（适配论文/报告）
end

%% Lambda
close all;

f1 = trajPlot.lambda;
t = trajPlot.t;
name = {'Horizontal ground contact force', 'Vertical ground contact force'};

% Data Phase 1 
time_1 = t;

% Data Phase 2
end_t = time_1(end);
time_2 = end_t + time_1;

f2 = f1;

for i = 1:2
    % Figure config
    figure(i);
    hold on;grid on;
    set(gcf, 'Position', [100, 100, 800, 300]);

    data1 = f1(:,i);
    data2 = f2(:,i);
    plot(time_1, data1, 'LineWidth', 2, 'DisplayName', 'On left leg')
    plot(time_2, data2, 'LineWidth', 2, 'DisplayName', 'On right leg')

    xlabel('Time [s]');
    ylabel('Force [N]');

    title(name{i});
    legend('show'); % 显示图例以区分颜色
    
    savefig(figure(i), sprintf('%s_v%g_N%g.fig', name{i},v,N));% save .fig

    filename_png = sprintf('%s_v%g_N%g.png', name{i},v,N);% save .png
    exportgraphics(gcf, filename_png, ...
    'Resolution', 300, ...       % 300DPI高清
    'BackgroundColor', 'white'); % 白色背景（适配论文/报告）
end

%% Cost u_squard
u = trajPlot.u;
cost = sum(sum(u.^2));

%% v2k_map
load summary.mat
%%
v = summary(2,:);
k_h = summary(3,:);
k_k = summary(4,:);

% Figure config
figure;
hold on;grid on;
% set(gcf, 'Position', [100, 100, 800, 300]);

plot(v, k_h, 'LineWidth', 2, 'DisplayName', 'k_h')
plot(v, k_k, 'LineWidth', 2, 'DisplayName', 'k_k')

xlabel('Average Velocity [m/s]');
ylabel('Stiffness [N·m/rad]');
% xlim([0.1,1.1]);
ylim([-1,21]);
title('Variation of Optimized Stiffness with Average Velocity');
legend('show'); % 显示图例以区分颜色

%%
% traj1 = T.x{1};
% traj2 = T.x{2};
%
% % 2. 创建图形窗口
% figure;
% hold on;
% grid on;
% 
% % 3. 绘制第一行轨迹 (红色，实线，加粗)
% plot(traj1(:, 1), traj1(:, 2), 'r-', 'LineWidth', 2, 'DisplayName', '轨迹 1 (Row 1)');
% 
% % 4. 绘制第二行轨迹 (蓝色，虚线，加粗)
% plot(traj2(:, 1), traj2(:, 2), 'b--', 'LineWidth', 2, 'DisplayName', '轨迹 2 (Row 2)');
% 
% % 5. 修饰图表
% xlabel('Time [s]');
% ylabel('State');
% title('轨迹曲线对比图');
% legend('show'); % 显示图例以区分颜色
% 
% %% Generate 
% x = [];
% u = [];
% t = [];
% p = [];
% lambda = [];
% 
% for i = 1:size(traj_full)
%     x = [x;traj_full(i).x];
%     u = [u;traj_full(i).u];
%     t = [t;traj_full(i).t];
%     p = [p;traj_full(i).p];
%     lambda = [lambda;traj_full(i).lambda];
% end
% traj.x = x;
% traj.u = u;
% traj.t = t;
% traj.p = p;
% traj.lambda = lambda;
% %%
% plot_traj(config, traj, [], [], []);
% 
% %%
% t = linspace(0, 10, 1000); 
% alpha = 0.2 + 0.1*sin(2*pi*0.5*t);
% 
% figure('Color','w','Position',[100,100,800,600]); % 创建白色背景的图窗，设置大小
% hold on; grid on; % 保留图层、显示网格（新手必加，提升可读性）
% 
% plot(t, alpha, 'r-', 'LineWidth',1.5, 'DisplayName','髋角α (rad)');
% 
% xlabel('时间 t (s)','FontSize',12); % x
% ylabel('数值','FontSize',12); % y
% title(name,'FontSize',14,'FontWeight','bold'); % title
% legend('Location','best','FontSize',10); % label
% xlim([0,10]); % x-limit
% 
% hold off;
% 
% %% 
