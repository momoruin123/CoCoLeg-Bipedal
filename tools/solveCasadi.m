function [solution, stats] = solveCasadi(config, Z, P, r, h, w, bounds, traj)
% SOLVECASADI Solve nonlinear optimization problem using CasADi and IPOPT
%
%   Forms and solves a nonlinear programming problem with dynamics constraints,
%   inequality constraints, and cost function using the specified trajectory
%   for initialization.
%
%   Inputs:
%     config  - system configuration structure
%     Z       - symbolic decision variables
%     P       - symbolic model parameters
%     r       - equality constraint residuals (dynamics + transitions)
%     h       - inequality constraint expressions
%     w       - cost function to minimize
%     bounds  - inequality constraint bounds [lower, upper]
%     traj    - trajectory data for initialization (struct or vector)
%
%   Outputs:
%     solution - solver solution structure containing optimal variables
%     stats    - solver statistics and performance metrics
%
%   author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

import casadi.*

% Initialize decision variables from trajectory data
if isstruct(traj)
    % Interpolate trajectory to get initial guess for optimization variables
    Z0 = interpolate_traj2Z(config, traj, config.N);
else
    % Use provided vector directly as initial guess
    Z0 = traj;
end

% Get numerical values for model parameters
P0 = getModelParameters(config, 0, 1);

% Formulate nonlinear programming problem
g = [r; h];  % Combine equality and inequality constraints
nlp = struct('x', Z, 'p', P, 'f', w, 'g', g);

opts = struct;

opts.ipopt.print_level = config.optConfig.print_level;
opts.ipopt.print_user_options = 'no';
opts.ipopt.sb = 'yes';
opts.print_time = false;

opts.ipopt.linear_solver = 'mumps';
opts.ipopt.nlp_scaling_method = 'gradient-based';

switch config.optConfig.tolerance
    case 1
        opts.ipopt.tol = 1e-5;
        opts.ipopt.jacobian_approximation = 'exact';
        opts.ipopt.hessian_approximation = 'limited-memory';
        opts.ipopt.max_iter = 1500;
    case 2
        opts.ipopt.tol = 1e-7;
        opts.ipopt.jacobian_approximation = 'exact';
        opts.ipopt.hessian_approximation = 'limited-memory';
        opts.ipopt.max_iter = 1500;
end

%%
% Create NLP solver instance
S = nlpsol('S', 'ipopt', nlp, opts);

% Set constraint bounds
% Equality constraints: r = 0, Inequality constraints: lb <= h <= ub
if isempty(bounds)
    lowerBound = [zeros(size(r))];  % [0; lb] for [r; h]
    upperBound = [zeros(size(r))];  % [0; ub] for [r; h]
else
    lowerBound = [zeros(size(r)); bounds(:,1)];  % [0; lb] for [r; h]
    upperBound = [zeros(size(r)); bounds(:,2)];  % [0; ub] for [r; h]
end

% Set States bounds
[lbx, ubx] = buildVariableBounds(config, Z);
% Solve the optimization problem
% solution = S('x0', Z0, 'lbg', lowerBound, 'ubg', upperBound, 'p', P0);

solution = S('x0', Z0, 'lbx', lbx ,'ubx', ubx, 'lbg', lowerBound, 'ubg', upperBound, 'p', P0);

% Extract solver statistics
stats = S.stats();
end

%% 
function [lbZ, ubZ] = buildVariableBounds(config, Z)
% Get informations from config
nPhase  = numel(config.phaseSequence);
na      = numel(config.optParameterNames) + config.optimizeTimeFlag;
N       = config.N;
nZ      = numel(Z);

% Optimization bounds of Z
q_lb  = config.q_lb;
q_ub  = config.q_ub;
dq_lb = config.dq_lb;
dq_ub = config.dq_ub;

input_lb    = config.input_lb;
input_ub    = config.input_ub;
dt_lb       = config.dt_bounds(1);
dt_ub       = config.dt_bounds(2);

% Optimization bounds of params
if numel(config.optParameterNames) > 0
    optParams_lb = config.optParameterLowerBound;                   % Minimum leg stiffness
    optParams_ub = config.optParameterUpperBound;
end

% Initialization
get_sizes   = cell(1, nPhase);
n_rows      = zeros(1, nPhase);
lbZ         = [];          % upperbounds
ubZ         = [];          % lowerbounds

for i = 1:nPhase
    % Get 'get_sizes()' function handle to each phase
    get_sizes_i = str2func([config.model_name, '.', config.phaseSequence{i}, '.get_sizes']);
    [nX, nq, nu, ~] = get_sizes_i();

    % Get the number of rows of Z corresponding to each phase
    n_rows_i = (N(i)+1) * (nX+nu+na);
    n_optValues = nX+nu+na;
    
    % Initialization
    lbZ_i = -inf(n_rows_i, 1);          % upperbounds
    ubZ_i =  inf(n_rows_i, 1);          % lowerbounds

    % Extend boundsX
    nBounds = numel(q_lb);
    if nq > nBounds
        q_lb  = [-inf*ones(nq-nBounds,1);q_lb];
        q_ub  = [ inf*ones(nq-nBounds,1);q_ub];
        dq_lb = [-inf*ones(nq-nBounds,1);dq_lb];
        dq_ub = [ inf*ones(nq-nBounds,1);dq_ub];
    end

    % Loop through the bounds of every collocation point
    for k = 1:N(i)+1
        offset  = (k-1) * n_optValues;

        % State bounds
        lbZ_i(offset + (1:nq)) = q_lb;
        ubZ_i(offset + (1:nq)) = q_ub;
        lbZ_i(offset + nq + (1:nq)) = dq_lb;
        ubZ_i(offset + nq + (1:nq)) = dq_ub;

        % Input bounds
        lbZ_i(offset + nX + (1:nu)) = input_lb;
        ubZ_i(offset + nX + (1:nu)) = input_ub;

        % Time step bound
        lbZ_i(offset + nX + nu + 1) = dt_lb;
        ubZ_i(offset + nX + nu + 1) = dt_ub;

        % Optimization params bound
        if numel(config.optParameterNames) > 0
            lbZ_i(offset + nX + nu + (2:na)) = optParams_lb;
            ubZ_i(offset + nX + nu + (2:na)) = optParams_ub;
        end

    end
    lbZ = [lbZ; lbZ_i];
    ubZ = [ubZ; ubZ_i];

end

end
