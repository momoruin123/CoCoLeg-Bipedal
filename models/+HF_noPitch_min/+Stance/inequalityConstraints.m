function [h, lb, ub, boundName] = inequaliyConstraints(config, Z, P_sym, N)
% INEQUALITYCONSTRAINTS Define inequality constraints for stance phase optimization
%
%   Inputs:
%     config - system configuration structure
%     Z      - decision variables vector
%     P_sym  - symbolic model parameters
%     N      - number of collocation points
%
%   Outputs:
%     h        - inequality constraint expressions in CASADI Symblic format
%     lb, ub   - lower and upper bounds for constraints
%     boundName- names of each constraint element
%
%   Author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*
    import HF_noPitch_min.Stance.*;  % Import stance phase functions
    
    % Load event function for contact force calculation
    [eventFcn, ~, ~] = S2F.createEventJump();     

    % Determine number of auxiliary variables
    na = numel(config.optParameterNames) + config.optimizeTimeFlag;

    % Get model dimensions
    [nx, ~, nu, ~] = get_sizes();

    % Initialize additional constraint counters
    n_add = 0;
    add_lb = [];
    add_ub = [];
    add_names = {};
    
    % Reshape decision variables into matrix form
    n_rows = nx + nu + na;
    n_columns = numel(Z) / n_rows;
    assert(mod(n_columns, 1) == 0 && n_columns > 0, 'n_columns must be a positive integer.');

    % Load numerical parameter values
    P_num = getModelParameters(config, 0, 1);
    
    % Initialize parameter bound structures
    paramValues = struct();
    param_lb = struct();
    param_ub = struct();

    % Get all parameter names from model
    getParameterNames = str2func([config.model_name, '.getParameterNames']);
    allParamNames = getParameterNames(config);
    
    % Set parameter values and bounds
    for i = 1:numel(allParamNames)
        paramName = allParamNames{i};
        isOptimized = ismember(paramName, config.optParameterNames);
        
        % Check if parameter has inequality bounds override
        hasIneqBound = isfield(config, 'inequalityParamBoundName') && ...
                       ismember(paramName, config.inequalityParamBoundName);
        
        if isOptimized
            % Use optimization bounds for optimized parameters
            idxOpt = find(strcmp(paramName, config.optParameterNames));
            
            if hasIneqBound
                % Override with inequality bounds
                idxIneq = find(strcmp(paramName, config.inequalityParamBoundName));
                param_lb.(paramName) = config.inequalityParamLowerBound(idxIneq);
                param_ub.(paramName) = config.inequalityParamUpperBound(idxIneq);
            else
                % Use standard optimization bounds
                param_lb.(paramName) = config.optParameterLowerBound(idxOpt);
                param_ub.(paramName) = config.optParameterUpperBound(idxOpt);
            end
            
            paramValues.(paramName) = config.optParameterLowerBound(idxOpt);
        else
            % Use fixed values for non-optimized parameters
            idx = getParameterIndex(paramName, config);
            paramValues.(paramName) = P_num(idx);
            
            if hasIneqBound
                % Apply inequality bounds even to fixed parameters
                idxIneq = find(strcmp(paramName, config.inequalityParamBoundName));
                param_lb.(paramName) = config.inequalityParamLowerBound(idxIneq);
                param_ub.(paramName) = config.inequalityParamUpperBound(idxIneq);
            else
                % Fixed parameters have equal lower/upper bounds
                param_lb.(paramName) = P_num(idx);
                param_ub.(paramName) = P_num(idx);
            end
        end
    end

    % Calculate natural frequency bounds for time step constraints
    eigenfreq_min = 1/(2*pi) * sqrt(param_lb.k_l / paramValues.m);  % [Hz]
    eigenfreq_max = 1/(2*pi) * sqrt(param_ub.k_l / paramValues.m);  % [Hz]
    T_max = 1 / eigenfreq_min; % Maximum period [s]
    T_min = 1 / eigenfreq_max; % Minimum period [s]

    % Initialize constraint arrays
    h = MX.zeros(2 * numel(Z));
    idx_h = 0;

    %% State Bounds
    % Hip angle bounds
    x_lb.alpha = -pi/2;
    x_ub.alpha = pi/2;
    
    % Leg length bounds (relative to natural length)
    x_lb.l = 1/2 * paramValues.l_0;
    x_ub.l = 3/2 * paramValues.l_0;

    % Angular velocity bounds 
    x_lb.dalpha = -2 * pi / T_min;  % [rad/s]
    x_ub.dalpha = 2 * pi / T_min;   % [rad/s]
   
    % Leg extension velocity bounds 
    x_lb.dl = -sqrt(5 * paramValues.g * x_ub.l);   
    x_ub.dl = sqrt(5 * paramValues.g * x_ub.l);

    %% Input Bounds
    epsilon_h = inf;  % Hip torque scaling factor
    epsilon_l = inf;  % Leg force scaling factor
  
    % Hip torque bounds (scaled by gravitational torque)
    u_lb.u_alpha = -epsilon_h * (1/2 * paramValues.m * paramValues.g * x_ub.l);
    u_ub.u_alpha = epsilon_h * (1/2 * paramValues.m * paramValues.g * x_ub.l);  

    % Leg force bounds (scaled by body weight)
    u_lb.u_l = -epsilon_l * paramValues.m * paramValues.g;  
    u_ub.u_l = epsilon_l * paramValues.m * paramValues.g;
    
    %% Auxiliary Variable Bounds
    % Time step bounds
    epsilon = 1e-2 / N;  % Minimum time step [s]
    a_lb.dt = epsilon;
    a_ub.dt = 3 * T_max / N;  % Maximum time step

    % Optimization parameter bounds
    for i = 1:numel(config.optParameterNames)
        paramName = config.optParameterNames{i};
        a_lb.(paramName) = config.optParameterLowerBound(i);
        a_ub.(paramName) = config.optParameterUpperBound(i);
    end

    %% Positive Contact Force Constraint
    % Extract variables from decision vector
    Z = reshape(Z, n_rows, n_columns);
    x = Z(1:nx, :);
    u = Z(nx+1:nx+nu, :);
    a = Z(nx+nu+1:n_rows, :);
    
    % Construct full parameter vector
    p_matrix = repmat(P_sym, 1, size(a,2));
    p_aux = a(end-na+2:end, :);
    p_full = [p_matrix; p_aux];  
    
    % Calculate contact forces
    pos_lambda = eventFcn(x, u, p_full);

    % Enforce positive contact force
    add_lb.lambda = -1e-4;  % Small tolerance for numerical stability
    add_ub.lambda = inf;
    
    % Extend variable set with contact forces
    Z_extended = [Z; pos_lambda];
    n_add = n_add + 1; 

    %% Jerk Constraints
    jerk_bound = inf;  % Jerk limit
    
    % Calculate control input derivatives (jerk)
    u_diff = [MX.zeros(2,1), diff(u')'];
    dt = repmat(a(1,:), size(u,1), 1);
    jerk = u_diff ./ dt;

    % Set jerk bounds for both controls
    add_lb.jerk_alpha = -jerk_bound;
    add_ub.jerk_alpha = jerk_bound;
    add_lb.jerk_l = -jerk_bound;
    add_ub.jerk_l = jerk_bound;

    % Extend variable set with jerk terms
    Z_extended = [Z_extended; jerk];
    n_add = n_add + 2; 

    %% putting it all together
    n_Z_extended = n_columns * (n_rows + n_add);
    Z_extended = reshape(Z_extended, n_Z_extended, 1);

    % Convert bound structures to vectors
    [x_lb, x_lb_names] = s2v(x_lb);
    [x_ub, ~] = s2v(x_ub);
    [u_lb, u_lb_names] = s2v(u_lb);
    [u_ub, ~] = s2v(u_ub);
    [a_lb, a_lb_names] = s2v(a_lb);
    [a_ub, ~] = s2v(a_ub);
    [add_lb, add_lb_names] = s2v(add_lb);
    [add_ub, ~] = s2v(add_ub);

    % Create bound vectors for all collocation points
    lb = repmat([x_lb; u_lb; a_lb(1:na); add_lb], n_columns, 1);
    ub = repmat([x_ub; u_ub; a_ub(1:na); add_ub], n_columns, 1);
    boundName = repmat([x_lb_names'; u_lb_names'; a_lb_names(1:na)'; add_lb_names'], n_columns, 1);

    % Remove unbounded constraints for efficiency
    unnecessaryConst = (lb == -inf) & (ub == inf);
    casadiIdx = find(unnecessaryConst == 0);
    lb = lb(~unnecessaryConst);
    ub = ub(~unnecessaryConst);
    boundName = boundName(~unnecessaryConst);

    % Assign constraint expressions
    h(idx_h+1:numel(lb)) = Z_extended(casadiIdx);
    idx_h = idx_h + numel(lb);

    % Truncate to actual number of constraints
    h = h(1:idx_h);
    boundName = boundName(1:idx_h);
end