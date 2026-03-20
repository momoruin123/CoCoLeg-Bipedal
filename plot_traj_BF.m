function plot_traj_BF(config, trajPlot)
%% config
nx = 10;
nq = 5;
nu = 4;

v = 0.7;
N = 150;

%% States
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

% Figure config
figure(1);
hold on;grid on;
set(gcf, 'Position', [100, 100, 800, 300]);

for i = 1:nx
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
    
    filename_png = sprintf('%s_v%g_N%g.png', name{i},v,N);% save .png
    exportgraphics(gcf, filename_png, ...
    'Resolution', 300, ...       % 300DPI高清
    'BackgroundColor', 'white'); % 白色背景（适配论文/报告）
end

%% Input
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

% Figure config
figure(2);
hold on;grid on;
set(gcf, 'Position', [100, 100, 800, 300]);

for i = 1:nu
    data1 = u1(:,i);
    data2 = u2(:,i);
    plot(time_1, data1, 'LineWidth', 2, 'DisplayName', 'On left leg')
    plot(time_2, data2, 'LineWidth', 2, 'DisplayName', 'On right leg')

    xlabel('Time [s]');
    ylabel('Torque [N·m]');

    title(name{i});
    legend('show'); % 显示图例以区分颜色
    
    filename_png = sprintf('%s_v%g_N%g.png', name{i},v,N);% save .png
    exportgraphics(gcf, filename_png, ...
    'Resolution', 300, ...       % 300DPI高清
    'BackgroundColor', 'white'); % 白色背景（适配论文/报告）
end

%% Lambda
f1 = trajPlot.lambda;
t = trajPlot.t;
name = {'Horizontal ground contact force', 'Vertical ground contact force'};

% Data Phase 1 
time_1 = t;

% Data Phase 2
end_t = time_1(end);
time_2 = end_t + time_1;

f2 = f1;

% Figure config
figure(3);
hold on;grid on;
set(gcf, 'Position', [100, 100, 800, 300]);

for i = 1:2
    data1 = f1(:,i);
    data2 = f2(:,i);
    plot(time_1, data1, 'LineWidth', 2, 'DisplayName', 'On left leg')
    plot(time_2, data2, 'LineWidth', 2, 'DisplayName', 'On right leg')

    xlabel('Time [s]');
    ylabel('Force [N]');

    title(name{i});
    legend('show'); % 显示图例以区分颜色
    
    filename_png = sprintf('%s_v%g_N%g.png', name{i},v,N);% save .png
    exportgraphics(gcf, filename_png, ...
    'Resolution', 300, ...       % 300DPI高清
    'BackgroundColor', 'white'); % 白色背景（适配论文/报告）
end
