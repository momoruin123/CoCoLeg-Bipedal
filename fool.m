clear;clc;
import BF_m.SingleStance.*
import casadi.*

%% Configeration of scrip
configfile_BF_m;

%% Prepration
L = 0.5;    % Steplength
T = 0.7;      % Steptime
N = 50;     % Segments
hk = T/N;   % Time interval between each Segment
ub_state = [0;  inf;  inf; -0.01; -0.01];
lb_state = [0; -inf; -inf; -inf; -inf];
% ub_state = [0;  0.5*pi;  0.5*pi ;      0;       0];
% lb_state = [0; -0.5*pi; -0.5*pi; -0.6*pi; -0.6*pi];
% ub_state = [0;  0.5*pi;  0.5*pi;  -0.1*pi;  -0.1*pi];
b_u      = [inf; inf; inf; inf];
% b_u      = [100; 100; 50; 50];

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
% linear interpolation
q_init = zeros(nq, N+1);
for k = 1:N+1
    t_k = (k-1)*hk;
    q_init(:,k) = q0_guess + (t_k/T)*(qT_guess - q0_guess);
end

% Initial guess of dq
% Average violocity
dq_init = (qT_guess - q0_guess)/T * ones(1, N+1);

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
for k = 1:N
    % State of k
    index_q = idx_q(k);
    index_u = idx_u(k);
    for i = 1:nq
        lb_num(index_q(i)) = lb_state(i);
        ub_num(index_q(i)) = ub_state(i);      
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
x0 = [q0; dq0];  
xT = [qT; dqT];
[~, jumpMap,~] = S2S.createEventJump();
x0_expected = jumpMap(xT, p);          
constraints_eq = vertcat(constraints_eq, x0 - x0_expected);

% % Pitch constraints
% for k = 1:N
%     % State of k
%     q_k = opt_vars(idx_q(k));
%     dq_k = opt_vars(idx_dq(k));
%     constraints_eq = vertcat(constraints_eq, q_k(1) - 0, dq_k(1) - 0);
% end

% Step length constraints
SwingFoot = cons_SwingFoot(qT, p);
constraints_eq = vertcat(constraints_eq, SwingFoot(1) - L, SwingFoot(2) - 0);

% --------------------------Inequation--------------------------
% Swing foot constraints
% Lift off: move upwards
Jacobian = J_SwingConstraint(x0, p);
dSwingFootPt = Jacobian * dq0;
constraints_ineq = vertcat(constraints_ineq, dSwingFootPt(2) - 0);
% Touch down: move downwards
Jacobian = J_SwingConstraint(xT, p);
dSwingFootPt = Jacobian * dqT;
constraints_ineq = vertcat(constraints_ineq, -dSwingFootPt(2) - 0);

% Swing foot >=0
for k = 1:N+1
    q_k = opt_vars(idx_q(k));
    SwingFoot = cons_SwingFoot(q_k, p);
    constraints_ineq = vertcat(constraints_ineq, SwingFoot(2) - 0);
end

% % Lambda constraints
% for k = 1:N+1
%     q_k = opt_vars(idx_q(k));
%     dq_k = opt_vars(idx_dq(k));
%     state = [q_k;dq_k];
%     [~, lambda] = BF_floating.left_Stance.dynamics(state, p);
%     constraints_ineq = vertcat(constraints_ineq, lambda(2) - 0);
% end

constraints = vertcat(constraints_eq, constraints_ineq);
lbg = [zeros(length(constraints_eq), 1); zeros(length(constraints_ineq), 1)];
ubg = [zeros(length(constraints_eq), 1); inf(length(constraints_ineq), 1)];

%% Solve
nlp = struct('x', opt_vars, 'f', J, 'g', constraints);

opts = struct;

% Output control
% Option 1:
opts.ipopt.print_level = 5;
opts.ipopt.tol = 1e-8;
opts.ipopt.dual_inf_tol = 1e-8;      
opts.ipopt.constr_viol_tol = 1e-6;
opts.ipopt.compl_inf_tol = 1e-8;     
opts.ipopt.linear_solver = 'mumps';  
opts.ipopt.jacobian_approximation = 'exact';  
opts.ipopt.hessian_approximation = 'exact';   
opts.ipopt.line_search_method = 'filter';
opts.ipopt.accept_every_trial_step = 'yes';  % more trial step
opts.ipopt.max_resto_iter = 10;
opts.ipopt.max_iter = 1000;
solver = nlpsol('solver', 'ipopt', nlp, opts);

% Option 2:
% opts.ipopt.print_level = 5;    
% opts.ipopt.tol = 1e-8;         
% opts.ipopt.constr_viol_tol = 1e-8;
% opts.ipopt.max_iter = 1000;    
% solver = nlpsol('solver', 'ipopt', nlp, opts);

sol = solver('x0', x0_num, ...
             'lbx', lb_num, ...
             'ubx', ub_num, ...
             'lbg', lbg, ...
             'ubg', ubg);

%% Optimization results
opt_vars_val = full(sol.x);

q_opt = reshape(opt_vars_val(1:nq*(N+1)), nq, N+1);      
dq_opt = reshape(opt_vars_val(nq*(N+1)+1:nX*(N+1)), nq, N+1); 
u_opt = reshape(opt_vars_val(nX*(N+1)+1:end), nu, N+1);     
t_grid = linspace(0, T, N+1);
x = [q_opt',dq_opt'];
traj.x = x;

%% Animation
save_flag = true;
file_name = mfilename;
time = datetime('now', 'Format','uuuuMMdd_HHmmss');
time = char(time);
gif_name = [file_name,'_',time];

getAnimationRABBIT(traj, 0, config, gif_name);


%% Analysis
import BF_floating.left_Stance.*

for i = 1:size(x,1)
    state = x(i,:);
    q  = state(1:nq);
    dq = state(nq+1:end);
    q  = [0,0,q];
    dq = [0,0,dq];
    state_f = [q, dq]';
    [~, lambda] = BF_floating.left_Stance.dynamics(state_f,p);
    if lambda < 0
        print('lambda was negative!!!')
        break;
    end
end
