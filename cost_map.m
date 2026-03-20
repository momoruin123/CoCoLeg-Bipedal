% run a cost map
clear;
clc;
%% BF_min: Hopping in place with minimum coordiantes
configfile_BF_min;
%%
load trajPlot_20260120_232037.mat
traj_init = trajPlot;

%% V_avg map
cost = [];
k_h = [];
init_v = 0.1;
% grid  = 0.01;
N = config.N;
for i = 1:91
    v_avg = init_v + (i-1)*grid;
    k_h = [k_h, v_avg];
    config.operatingCond.v_avg = v_avg;
    [Z, P, w, r, constraintNames, h, bounds, boundNames] = getConstraintsAndCosts(config, traj_init);
    [solutionIPOPT,  stats] = solveCasadi(config, Z, P, r, h, w, bounds, traj_init);

    % Unpack result
    solution = full(solutionIPOPT.x);
    u_idx = 11;
    u = [];
    for j = 1:N+1
        u_temp = solution((j-1)*15+(11:14));
        u = [u,u_temp];
    end
    u = u';
    cost = [cost, sum(sum(u.^2))/26];
end


%%
load cost_v_map.mat

figure;
hold on; 
grid on;
% set(gcf, 'Position', [100, 100, 800, 300]);

cost = cost_v_map(1,:);
k_h = cost_v_map(2,:);

plot(k_h, cost, 'LineWidth', 2);

title('Variation of Squared Cost with Average Velocity');
xlabel('average vel [m/s]');
ylabel('cost');

%% k map
cost = [];
k_h = [];
k_k = [];

init_k_h = 7.10633207498287;
init_k_k = 0.663507536778245;

grid_k_h = 1;
grid_k_k = 0.1;


N = config.N;
numOptValue = 15;
u_idx = 11;

for i = 1:5
    for k = 1:5
        disp(i)

        config.k_h = init_k_h + (i-1)*grid_k_h;
        config.k_k = init_k_k + (k-1)*grid_k_k;
        k_h(i,k) = config.k_h;
        k_k(i,k) = config.k_k;

        [Z, P, w, r, constraintNames, h, bounds, boundNames] = getConstraintsAndCosts(config, traj_init);
        [solutionIPOPT,  stats] = solveCasadi(config, Z, P, r, h, w, bounds, traj_init);

        % Unpack result
        solution = full(solutionIPOPT.x);
        f = solutionIPOPT.f;
        cost(i,k) = full(f);
    end
end
%%
save('cost_k_5.mat',"cost")
save('k_h_5.mat',"k_h")
save('k_k_5.mat',"k_k")

%% 加载数据（使用时取消注释）
% load summary_k.mat
a1 = k_h;    % 6×6矩阵 (k_h)
a2 = k_k;    % 6×6矩阵 (k_k)
a3 = cost;   % 6×6矩阵 (cost)

% 绘制三维曲面图（推荐用surf，视觉效果更优）
figure('Color','w');
surf(a1, a2, a3);  % 直接用6×6矩阵绘制曲面，无需转一维
% 可选：如果想要网格线+填充色的组合效果，用surf；如果只要网格线，用mesh(a1,a2,a3)

% 图形美化与标注
xlabel('k_h'); 
ylabel('k_k'); 
zlabel('Cost');
title('Variation of Squared Cost with Stiffness');
colormap(jet);     % 配色方案，可替换为parula/gray/hot等
colorbar;          % 显示颜色条（对应cost值）
shading interp;    % 平滑着色，消除网格锯齿，让曲面更细腻
grid on;           % 显示网格
rotate3d on;       % 允许鼠标旋转视图
alpha(0.9);        % 曲面透明度（0-1，可选）

%% 生成随机三维数据
x = randn(1000,1);  % 1000个随机x
y = randn(1000,1);  % 1000个随机y
z = ones(1000,1);  % 1000个随机z

figure('Color','w');
scatter3(x, y, z, 10, z, 'filled');  % 10是点的大小，z表示按z值配色，filled填充点
xlabel('X轴'); ylabel('Y轴'); zlabel('Z轴');
title('三维散点图 (scatter3)');
colormap(jet);
colorbar;
grid on;