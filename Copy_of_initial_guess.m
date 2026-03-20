function traj_guess = initial_guess(config, init_stepTime)
traj_guess = struct();
nq = 5;
nu = 4;
N  = config.N;           % Extract segment number
hk = init_stepTime/N;   % Time interval between each Segment% Initial guess of q
init_x = config.x0;
% Transformation matrix
Trans_matrix = [1 0 0 0 0;
    0 0 1 0 0;
    0 1 0 0 0;
    0 0 0 0 1;
    0 0 0 1 0;];

% Initial guess of q
q0_guess = init_x(1:nq);
qT_guess = Trans_matrix * q0_guess;

dq0_guess = init_x(nq+1:end);
dqT_guess = Trans_matrix * dq0_guess;

% linear interpolation
t_init = zeros(1, N+1);
q_init = zeros(nq, N+1);
for k = 1:N+1
    t_k = (k-1)*hk;
    t_init(k) = t_k;
    q_init(:,k) = q0_guess + (t_k/init_stepTime)*(qT_guess - q0_guess);
end

dq_init = zeros(nq, N+1);
for k = 1:N+1
    t_k = (k-1) * hk;
    dq_init(:,k) = dq0_guess + (t_k/init_stepTime) * (dqT_guess - dq0_guess);
end

% Initial guess of u
u_init = zeros(nu, N+1);

traj_guess.x = [q_init; dq_init]';
traj_guess.u = u_init';
traj_guess.t = t_init';

% Options: depends on optimizing stiffness or not.
% traj.p = [20*ones(N+1,1), 20*ones(N+1,1)];
end