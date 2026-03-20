function [r, w, constraintName] = dynamics_residuals(config, packageName, Z, P, N)
% DYNAMICS_RESIDUALS Compute dynamics constraints and cost using Hermite-Simpson collocation
%
%   Inputs:
%     config      - system configuration structure
%     packageName - name of model package containing dynamics functions
%     Z           - decision variables matrix
%     P           - model parameters
%     N           - number of collocation intervals
%
%   Outputs:
%     r            - dynamics constraint residuals
%     w            - integrated cost
%     constraintName- names of constraint elements
%
    import casadi.*
    
    % Get model dimensions from the specified package
    get_sizes = str2func([packageName, '.get_sizes']);

    % Extract collocation and interpolation function names from config, use defaults if not specified
    if isfield(config, 'collocationScheme')
        collocation_func = str2func([config.collocationScheme, '_collocation']);
    else
        % standard if not specified in config file
        collocation_func = @hermiteSimpson_collocation;
    end
    
    if isfield(config, 'inputInterpolation') 
        inputInterp_func =  str2func([config.inputInterpolation, '_input']);
    else
        % standard if not specified in config file
        inputInterp_func = @piecewiseLinear_input;
    end
    
    % Extract configuration parameters
    cost = config.costName;  % Type of cost function to use
    na   = numel(config.optParameterNames) + config.optimizeTimeFlag;  % Number of auxiliary variables

    % Get model dimensions
    [nx, ~, nu, ~] = get_sizes();

    % Create mapped function for piecewise linear control inputs
    % Uses average of control inputs at interval boundaries
    r_u_Fun = inputInterp_func(nu, N);

    % Create mapped function for auxiliary variables (constant over intervals)
    % Ensures auxiliary variables remain constant across collocation points
    a_k     = MX.sym('a_k', na, 1);
    a_k1    = MX.sym('a_k1', na, 1);
    r_a_fun = Function('r_a', {a_k, a_k1}, {-a_k1 + a_k});
    r_a_fun = r_a_fun.map(N);  % Vectorize over N intervals

    % Initialize outputs
    r = MX.zeros(numel(Z));        % Residual vector for dynamics constraints
    w = MX.zeros(1);               % Total cost
    constraintName = cell(numel(Z),1);  % Constraint names
    idx_r = 0;                     % Index counter for residuals

    % Reshape decision variables into structured format
    % Z contains: states (nx), controls (nu), auxiliary variables (na) for N+1 points
    Z = reshape(Z, nx+nu+na, N+1);
    x = Z(1:nx, :);      % State variables
    u = Z(nx+1:nx+nu, :); % Control inputs  
    a = Z(nx+nu+1:nx+nu+na, :); % Auxiliary variables (time step + parameters)
    if config.optimizeTimeFlag
        dt = a(1,:);         % Time steps for each interval
    else
        hk = config.Time / N;         % Time steps for each interval
        dt = hk * ones(1, N+1);   
    end

    % Create Hermite-Simpson collocation functions
    [x_mid_fun, r_dx_fun] = collocation_func(nx, N);

    % Extract variables at interval boundaries
    x_k   = x(:, 1:N);    % States at start of intervals
    x_k1  = x(:, 2:N+1);  % States at end of intervals
    u_k   = u(:, 1:N);    % Controls at start of intervals
    u_k1  = u(:, 2:N+1);  % Controls at end of intervals
    a_k   = a(:, 1:N);    % Auxiliary variables at start of intervals
    a_k1  = a(:, 2:N+1);  % Auxiliary variables at end of intervals

    % Compute dynamics at all points
    f = get_dynamics(packageName, x, u, a, P, na);
    f_k   = f(:,1:N);     % Dynamics at start of intervals
    f_k1  = f(:,2:N+1);   % Dynamics at end of intervals

    % Compute midpoint values using Hermite-Simpson scheme
    x_mid = x_mid_fun(x_k, x_k1, f_k, f_k1, dt(:, 1:N));  % Midpoint states
    u_mid = r_u_Fun(u_k, u_k1);                           % Midpoint controls
    f_mid = get_dynamics(packageName, x_mid, u_mid, a_k, P, na);  % Midpoint dynamics

    % Compute dynamics residuals using Hermite-Simpson defect constraints
    r_dx = r_dx_fun(x_k, x_k1, f_k, f_mid, f_k1, dt(:, 1:N));  % State defect constraints
    r_a  = r_a_fun(a_k, a_k1);                                 % Auxiliary variable constraints

    % Combine and flatten all residuals
    r_total = [r_dx; r_a];                          % Combine state and auxiliary residuals
    r_total = reshape(r_total, (nx+na)*N, 1);       % Flatten to vector
    
    % Store residuals and constraint names
    r(idx_r + (1:numel(r_total))) = r_total;
    constraintName(idx_r + (1:numel(r_total))) = {'dynamics residuals'};
    idx_r = idx_r + numel(r_total);

    % Compute cost using Simpson's rule integration
    w_k   = compute_cost_(packageName, cost, x_k, u_k, dt, N);     % Cost at start points
    w_mid = compute_cost_(packageName, cost, x_mid, u_mid, dt, N); % Cost at midpoints  
    w_k1  = compute_cost_(packageName, cost, x_k1, u_k1, dt, N);  % Cost at end points

    % Simpson's rule integration: (h/6)*[f(a) + 4f((a+b)/2) + f(b)]
    w = w + 1/6*(w_k + 4*w_mid + w_k1); 

    % Truncate outputs to actual sizes
    r = r(1:idx_r);
    constraintName = constraintName(1:idx_r);
end

function r_u_Fun = piecewiseLinear_input(nu, N)
    % piecewiseLinear_input Create mapped function for piecewise linear control interpolation
    %
    %   Inputs:
    %     nu - number of control inputs
    %     N  - number of intervals to vectorize over
    %
    %   Output:
    %     r_u_Fun - mapped function that computes average control at midpoints

    import casadi.*
    
    % Define symbolic variables for control inputs at interval boundaries
    u_k  = MX.sym('u_k', nu, 1);
    u_k1 = MX.sym('u_k1', nu, 1);
    
    % Create function that computes average control (piecewise linear interpolation)
    r_u_Fun = Function('r_u_Fun', {u_k, u_k1}, {(1/2)*(u_k + u_k1)});
    
    % Vectorize function over N intervals
    r_u_Fun = r_u_Fun.map(N);
end

function [x_mid_fun, r_dx_fun] = hermiteSimpson_collocation(nx, N)
% HERMITESIMPSON_COLLOCATION Create Hermite-Simpson collocation functions
%
%   Creates mapped functions for computing midpoint states and dynamics
%   defect constraints using the Hermite-Simpson collocation scheme
%
%   Inputs:
%     nx - number of states
%     N  - number of intervals to vectorize over
%
%   Outputs:
%     x_mid_fun - function to compute midpoint states
%     r_dx_fun  - function to compute dynamics defect constraints

    import casadi.*
    
    % Define symbolic variables
    x_k     = MX.sym('x_k', nx, 1);    % State at start of interval
    x_k1    = MX.sym('x_k1', nx, 1);   % State at end of interval  
    dx_k    = MX.sym('dx_k', nx, 1);   % Dynamics at start of interval
    dx_mid  = MX.sym('dx_mid', nx, 1); % Dynamics at midpoint
    dx_k1   = MX.sym('dx_k1', nx, 1);  % Dynamics at end of interval
    dt      = MX.sym('dt', 1);         % Time step

    % Hermite-Simpson midpoint state calculation
    % x_mid = 0.5*(x_k + x_k1) + (dt/8)*(dx_k - dx_k1)
    x_mid_fun = Function('x_mid_fun', {x_k, x_k1, dx_k, dx_k1, dt}, ...
                       {0.5*(x_k + x_k1) + (dt/8) * (dx_k - dx_k1)});
    
    % Hermite-Simpson defect constraint  
    % r_dx = dt*dx_mid + 1.5*(x_k - x_k1) + (1/4)*dt*(dx_k + dx_k1)
    r_dx_fun = Function('r_dx_fun', {x_k, x_k1, dx_k, dx_mid, dx_k1, dt}, ...
                        {dt*dx_mid + 1.5*(x_k - x_k1) + (1/4)*dt*(dx_k + dx_k1)});

    % Vectorize functions over N intervals
    x_mid_fun = x_mid_fun.map(N);
    r_dx_fun = r_dx_fun.map(N);
end

function f = get_dynamics(packageName, x, u, a, p, na) 
% GET_DYNAMICS Evaluate system dynamics at multiple points
%
%   Inputs:
%     packageName - name of model package
%     x           - state variables
%     u           - control inputs  
%     a           - auxiliary variables
%     p           - model parameters
%     na          - number of auxiliary variables
%
%   Output:
%     f           - dynamics evaluated at all points

    import casadi.*

    % Create dynamics function from model package
    createDynamics = str2func([packageName, '.createDynamics']);
    dyn = createDynamics();
    dyn = dyn.map(size(x,2));  % Vectorize over all points
    
    % Construct full parameter vector combining fixed and optimized parameters
    % When na = 2: a = [dt; p_opt], so p_aux = a(2:end) = p_opt
    % When na = 1: a = [dt], so p_aux = [] (no optimized parameters)
    p_matrix = repmat(p, 1, size(a,2));      % Fixed parameters repeated for all points
    p_aux    = a(end-na+2:end, :);          % Optimized parameters (empty if na=1)
    p_full   = [p_matrix; p_aux];           % Combined parameter vector
    
    % Evaluate dynamics at all points
    f = dyn(x, u, p_full);
end

function w_k = compute_cost_(packageName, cost, x_k, u_k, dt, N)
% COMPUTE_COST_ Compute cost function at given points
%
%   Inputs:
%     packageName - name of model package  
%     cost        - cost function type
%     x_k         - state variables
%     u_k         - control inputs
%     a_k         - auxiliary variables (contains time steps)
%
%   Output:
%     w_k         - computed cost

    % Get time steps and control dimension
    dt = dt(1,1:N);  % Time steps for integration
    nu = size(u_k, 1);  % Number of controls

    % Select and compute cost based on specified type
    if strcmpi(cost, 'none')
        w_k = 0;  % No cost
    elseif strcmpi(cost, 'u_squared')
        % Quadratic control cost: sum(dt * u^2)
        w_k = sum(dot(repmat(dt, nu, 1).*u_k, u_k)); 
    elseif strcmpi(cost, 'weighted_u_squared')
        % Weighted quadratic control cost with different weights for each control
        weights = [1; 1/4];  % Different weights for different controls
        w_k = sum(dot(repmat(dt, nu, 1).*u_k, weights.*u_k)); 
    elseif strcmpi(cost, 'positive_mechanical_work')
        % Positive mechanical work cost from model-specific function
        posMechWorkFcn = str2func([packageName, '.posMechWork']);
        w_k = sum(dt.*posMechWorkFcn(x_k, u_k));  % Integrate using time steps
    end
end