clear;clc;
import BF_min.SingleStance.*
import casadi.*

% 更新日志：
% 终于他妈的做出来了，问题是啥呢？问题是我后来加上了phi的约束，让它等于0，这可能对
% 在实现碰撞的同时完成周期步态是很难的或许说是不可能的。有意思，可能可以从侧面说明
% 人类最节能的走路方式是一定要上半身的phi有一定晃动的

%% Configeration of scrip
configfile_BF_min;

%% Prepration
L = 0.3;    % Steplength
T = 0.7;      % Steptime
N = 25;     % Segments
hk = T/N;   % Time interval between each Segment

ub_state = [0.5*pi;  0.5*pi;  0.5*pi; -0.01; -0.01];
lb_state = [-0.5*pi; -0.5*pi; -0.5*pi; -inf; -inf];

dqLow = -10*ones(5,1);
dqUpp = 10*ones(5,1);

% ub_state = [0;  inf;  inf; -0.01; -0.01];
% lb_state = [0; -inf; -inf; -inf; -inf];

% ub_state = [0;  0.5*pi;  0.5*pi ;      0;       0];
% lb_state = [0; -0.5*pi; -0.5*pi; -0.6*pi; -0.6*pi];

% ub_state = [0;  0.5*pi;  0.5*pi;  -0.1*pi;  -0.1*pi];
% lb_state = [0; -0.5*pi; -0.5*pi; -0.6*pi; -0.6*pi];

b_u      = [inf; inf; inf; inf];
% b_u      = [150; 150; 150; 150];

% Number of system
[nX, nq, nu, np] = get_sizes();
nv = (nX + nu)*(N + 1);  % Number of full variables

% Parameters of model
[p, paraName] = getFullParameters_biped(config);

%% Initialization
% Initialize null optimized variables
opt_vars = MX.sym('opt_vars',nv,1);

% Initial guess of q
q0_guess = [0;0.7;-0.5;-1;-0.1];
qT_guess = [0;-0.5;0.7;-0.1;-1];

dq0_guess = (qT_guess-q0_guess)/T;
dqT_guess = dq0_guess;

% linear interpolation
q_init = zeros(nq, N+1);
for k = 1:N+1
    t_k = (k-1)*hk;
    q_init(:,k) = q0_guess + (t_k/T)*(qT_guess - q0_guess);
end

dq_init = zeros(nq, N+1);
for k = 1:N+1
    t_k = (k-1) * hk;
    dq_init(:,k) = dq0_guess + (t_k/T) * (dqT_guess - dq0_guess);
end

% Initial guess of u
u_init = zeros(nu, N+1);

% Initial values ​​of the optimized variable
x0_num = [q_init(:); dq_init(:); u_init(:)];

% Initialize null optimized variables limit
lb_num = -inf(nv, 1);
ub_num = inf(nv, 1);

% Index functions
idx_q  = @(k) (k-1)*nq+1 : k*nq;          
idx_dq = @(k) (k-1)*nq+1 + (N+1)*nq : k*nq + (N+1)*nq;     
idx_u  = @(k) (k-1)*nu+1 + 2*(N+1)*nq : k*nu + 2*(N+1)*nq; 

%% Optimization 
% Objective function
J = 0;
for k = 1:N
    u_k = opt_vars(idx_u(k));
    u_k1 = opt_vars(idx_u(k+1));
    
    J = J + 0.5 * hk * (u_k'*u_k + u_k1'*u_k1);
end

% Constraints
constraints_eq = MX.zeros(0,1);  
constraints_ineq = MX.zeros(0,1);

% --------------------------State limit-----------------------
for k = 1:N+1
    % State of k
    index_q = idx_q(k);
    index_dq = idx_dq(k);
    index_u = idx_u(k);
    for i = 1:nq
        lb_num(index_q(i)) = lb_state(i);
        ub_num(index_q(i)) = ub_state(i);      
    end
    for i = 1:nq
        lb_num(index_dq(i)) = dqLow(i);
        ub_num(index_dq(i)) = dqUpp(i);      
    end
    for i = 1:nu
        lb_num(index_u(i)) = -b_u(i);
        ub_num(index_u(i)) =  b_u(i); 
    end
end

% --------------------------Equation--------------------------
% Dynamic constraints
Dyn = createDynamics();
for k = 1:N
    % State of k
    q_k = opt_vars(idx_q(k));
    dq_k = opt_vars(idx_dq(k));
    x_k = [q_k; dq_k];
    u_k = opt_vars(idx_u(k));
    % State of k+1
    q_k1 = opt_vars(idx_q(k+1));
    dq_k1 = opt_vars(idx_dq(k+1));
    x_k1 = [q_k1; dq_k1];
    u_k1 = opt_vars(idx_u(k+1));
    % Dynamic 
    dx_k = Dyn(x_k, u_k, p);
    dx_k1 = Dyn(x_k1, u_k1, p);
    % Constraints
    eq_q = q_k1 - q_k - 0.5 * hk * (dx_k(1:nq) + dx_k1(1:nq));
    eq_dq = dq_k1 - dq_k - 0.5 * hk * (dx_k(nq+1:end) + dx_k1(nq+1:end));
    
    constraints_eq = vertcat(constraints_eq, eq_q, eq_dq);
end

% Periodic constraints
q0 = opt_vars(idx_q(1));
dq0 = opt_vars(idx_dq(1));
qT = opt_vars(idx_q(N+1));
dqT = opt_vars(idx_dq(N+1));
u0 = opt_vars(idx_u(1));
uT = opt_vars(idx_u(N+1));

x0 = [q0; dq0];  
xT = [qT; dqT];
[~, jumpMap,~] = S2S.createEventJump();
x0_expected = jumpMap(xT, uT, p);   
constraints_eq = vertcat(constraints_eq, x0 - x0_expected);

% Step length constraints
% SwingFoot = cons_SwingFoot_m(qT, p);
% constraints_eq = vertcat(constraints_eq, SwingFoot(1) - L, SwingFoot(2) - 0);
% constraints_eq = vertcat(constraints_eq, SwingFoot(2) - 0);

% velocity constraints
v_target = 0.7;
linkPts_0 = LinkPositions([0;0;q0], p);
linkPts_T = LinkPositions([0;0;qT], p);
hipPt_0 = linkPts_0(:,3);
hipPt_T = linkPts_T(:,3);

v_real = (hipPt_T(1) - hipPt_0(1))/T;
constraints_eq = vertcat(constraints_eq, v_real - v_target);

% --------------------------Inequation--------------------------
% Swing foot velocity constraints
% Lift off: move upwards
Jacobian = J_SwingConstraint_m(x0, p);
dSwingFootPt = Jacobian * dq0;
constraints_ineq = vertcat(constraints_ineq, dSwingFootPt(2) - 0);
% Touch down: move downwards
Jacobian = J_SwingConstraint_m(xT, p);
dSwingFootPt = Jacobian * dqT;
constraints_ineq = vertcat(constraints_ineq, -dSwingFootPt(2) - 0);

% Swing foot >=0
for k = 2:N
    q_k = opt_vars(idx_q(k));
    SwingFoot = cons_SwingFoot_m(q_k, p);
    constraints_ineq = vertcat(constraints_ineq, SwingFoot(2) - 0);
end

% Lift foot up limit
mid_k = floor(N/2) + 1;
q_k = opt_vars(idx_q(mid_k));
SwingFoot = cons_SwingFoot_m(q_k, p);
constraints_eq = vertcat(constraints_eq, SwingFoot(2) - 0.1);

% Lambda constraints
slack = 1e-6;
for k = 1:N+1
    q_k = opt_vars(idx_q(k));
    dq_k = opt_vars(idx_dq(k));  
    input = opt_vars(idx_u(k));
    lambda = compute_lambda(q_k, dq_k, input, p);

    constraints_ineq = vertcat(constraints_ineq, lambda(2) - slack);
end

constraints = vertcat(constraints_eq, constraints_ineq);
lbg = [zeros(length(constraints_eq), 1); zeros(length(constraints_ineq), 1)];
ubg = [zeros(length(constraints_eq), 1); inf(length(constraints_ineq), 1)];
% 
% constraints = constraints_eq;
% lbg = [zeros(length(constraints_eq), 1)];
% ubg = [zeros(length(constraints_eq), 1)];
%% Solve
nlp = struct('x', opt_vars, 'f', J, 'g', constraints);

opts = struct;

% Output control
% Option 1
% opts.ipopt.print_level = 5;
% opts.ipopt.tol = 1e-8;
% opts.ipopt.dual_inf_tol = 1e-8;      
% opts.ipopt.constr_viol_tol = 1e-6;
% opts.ipopt.compl_inf_tol = 1e-8;     
% opts.ipopt.linear_solver = 'mumps';  
% opts.ipopt.jacobian_approximation = 'exact';  
% opts.ipopt.hessian_approximation = 'exact';   
% opts.ipopt.line_search_method = 'filter';
% opts.ipopt.accept_every_trial_step = 'yes';  % more trial step
% opts.ipopt.max_resto_iter = 10;
% opts.ipopt.max_iter = 2000;
% solver = nlpsol('solver', 'ipopt', nlp, opts);

% Option 2
% opts.ipopt.print_level = 5;    
% opts.ipopt.tol = 1e-8;         
% opts.ipopt.constr_viol_tol = 1e-8;
% opts.ipopt.max_iter = 2000;
% opts.ipopt.acceptable_tol = 1e-6; % 容许的精度，防止最后几步死磕
% opts.ipopt.mu_strategy = 'adaptive'; % 自适应势垒参数，对非线性强的系统有效
% opts.ipopt.nlp_scaling_method = 'gradient-based'; % 自动缩放，平衡不同量级的数据solver = nlpsol('solver', 'ipopt', nlp, opts);
% solver = nlpsol('solver', 'ipopt', nlp, opts);

% Option 3
% 极简配置（确保IPOPT能启动）
% opts = struct();
% opts.ipopt.print_level = 5;          % 核心参数：日志输出级别（0-12）
% opts.ipopt.tol = 1e-8;               % 核心参数：主收敛精度
% opts.ipopt.constr_viol_tol = 1e-8;   % 核心参数：约束违反容忍度
% opts.ipopt.max_iter = 2000;          % 核心参数：最大迭代次数
% opts.ipopt.acceptable_tol = 1e-6;    % 核心参数：容许精度
% 
% % 初始化求解器（验证基础配置是否有效）
% solver = nlpsol('solver', 'ipopt', nlp, opts);
% 
% sol = solver('x0', x0_num, ...
%              'lbx', lb_num, ...
%              'ubx', ub_num, ...
%              'lbg', lbg, ...
%              'ubg', ubg);

%%
% Relaxed Precision
opts1 = struct();
opts1.ipopt.print_level = 5;          
opts1.ipopt.tol = 1e-3;               
opts1.ipopt.constr_viol_tol = 1e-3;   
opts1.ipopt.max_iter = 1000;          
opts1.ipopt.acceptable_tol = 1e-2;    

solver1 = nlpsol('solver1', 'ipopt', nlp, opts1);

sol1 = solver1('x0', x0_num, ...
               'lbx', lb_num, ...
               'ubx', ub_num, ...
               'lbg', lbg, ...
               'ubg', ubg);

% High precision
opts2 = struct();
opts2.ipopt.print_level = 5;          
opts2.ipopt.tol = 1e-8;               
opts2.ipopt.constr_viol_tol = 1e-8;  
opts2.ipopt.max_iter = 2000;          
opts2.ipopt.acceptable_tol = 1e-6;    

solver2 = nlpsol('solver2', 'ipopt', nlp, opts2);

% Using the results of the first optimization as an initial guess
sol2 = solver2('x0', full(sol1.x), ...
               'lbx', lb_num, ...
               'ubx', ub_num, ...
               'lbg', lbg, ...
               'ubg', ubg);

% 最终结果
opt_vars_val = full(sol2.x);

%% Optimization results
% opt_vars_val = full(sol.x);

q_opt = reshape(opt_vars_val(1:nq*(N+1)), nq, N+1);      
dq_opt = reshape(opt_vars_val(nq*(N+1)+1:nX*(N+1)), nq, N+1); 
u_opt = reshape(opt_vars_val(nX*(N+1)+1:end), nu, N+1);     
t_grid = linspace(0, T, N+1);
x = [q_opt',dq_opt'];

traj.x = x;
traj.u = u_opt';
traj.t = t_grid';
% traj.p = p;

%% Animation
save_flag = true;
file_name = mfilename;
time = datetime('now', 'Format','uuuuMMdd_HHmmss');
time = char(time);
if save_flag
    gif_name = [file_name,'_',time];
    getAnimationRABBIT(traj, 0, config, gif_name);
else
    getAnimationRABBIT(traj, 0, config);
end

%% Save result
mat_name = [gif_name,'.mat'];
save(mat_name, 'traj');

%% Analysis
for i = 1:size(x,1)
    state = x(i,:);
    q  = state(1:nq)';
    dq = state(nq+1:end)';
    input = u_opt(:,i);
    lambda = compute_lambda(q, dq, input, p);
    if lambda(2) < -1e-6
        fprintf("lambda was negative!!!")
        disp(lambda(2));
        break;
    else
        fprintf("lambda was positive:%f\n", lambda(2));
    end
end

%% Calculate the contact force between the ground and the stance foot
function lambda = compute_lambda(q, dq, u, p)
    import BF_min.SingleStance.*
    temp = [0; 0];
    q = vertcat(temp, q);
    dq = vertcat(temp, dq);
    nq = numel(q);
    X = [q;dq];

    M_f = massMatrix_f(X, p);  
    F_f= forces_f(X, u, p);
    
    J_Stance = J_StanceConstraint_f(X, p);
    dJ_Stance = dJ_StanceConstraint_f(X, p);
    
    % gamma = dJ* dq
    gamma = dJ_Stance * dq;
    
    A = [M_f    , -J_Stance';
        J_Stance, zeros(size(J_Stance,1))];
    B = [F_f;
        -gamma];
    
    sol = A \ B;
    lambda = sol(nq+1:end);
end