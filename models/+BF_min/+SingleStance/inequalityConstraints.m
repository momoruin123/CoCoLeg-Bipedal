function [h, lb, ub, boundName] = inequalityConstraints(config, Z, P_sym, N)
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

    import casadi.*
    import BF_min.SingleStance.*;  % Import stance phase functions 

    % Determine number of auxiliary variables
    na = numel(config.optParameterNames) + config.optimizeTimeFlag;

    % Get model dimensions
    [nx, nq, nu, ~] = get_sizes();

    % Reshape decision variables into matrix form
    n_rows = nx + nu + na;
    n_columns = numel(Z) / n_rows;
    assert(mod(n_columns, 1) == 0 && n_columns > 0, 'n_columns must be a positive integer.');

    % Initialize constraint arrays
    h = MX(0, 1);
    lb = [];
    ub = [];
    boundName = cell(0, 1);

    %% Additional Dynamic Constraints
    % Extract variables from decision vector
    Z = reshape(Z, n_rows, n_columns);
    x = Z(1:nx, :);
    u = Z(nx+1:nx+nu, :);
    a = Z(nx+nu+1:n_rows, :);
    q = x(1:nq, :);
    dq = x(nq+1:end, :);

    % Construct full parameter vector
    p_matrix = repmat(P_sym, 1, size(a,2));
    p_aux = a(end-na+2:end, :);
    p_full = [p_matrix; p_aux];  
    
    %% 1. Ground Clearance Constraint   
    % Swing foot >=0
    swingFoot_pos = pos_SwingFoot_m(x, p_full);
    swingFoot_pos_y = swingFoot_pos(2,:)';
    ineq = swingFoot_pos_y;

    % Global clearance constraint
    global_lower_bound = 0;
    global_upper_bound = inf;
    name_bound = 'ground_clearance';

    % Bound
    lower_bound = global_lower_bound * ones(n_columns, 1);
    upper_bound = global_upper_bound * ones(n_columns, 1);

    % Local clearance constraint
    local_lower_bound = 0.8*0.01;
    offset = floor(n_columns*0.1);
    for k = offset:(n_columns-offset)
        lower_bound(k) = local_lower_bound;
    end  

    % Extend
    h = [h; ineq];
    lb = [lb; lower_bound];
    ub = [ub; upper_bound];
    for i = 1:n_columns
        boundName{end+1,1} = name_bound;
    end

    %% 2. Swing foot velocity constraints
    % Extract necessary state
    x0 = x(:, 1);
    dq0 = dq(:, 1);
    xT = x(:, end);
    dqT = dq(:, end);
    p_full_0 = p_full(:, 1);
    p_full_T = p_full(:, end);
    %% Lift off: Lift vertical velocity > 0
    Jacobian_lo = J_SwingConstraint_m(x0, p_full_0);
    v0_swingFoot = Jacobian_lo * dq0;
    ineq = v0_swingFoot(2);

    % Bounds
    lower_bound = 0;
    upper_bound = inf;
    name_bound = 'lift_off_vel_y';

    % Extend
    h = [h; ineq];
    lb = [lb; lower_bound];
    ub = [ub; upper_bound];
    boundName{end+1,1} = name_bound;
    
    %% Touch down: Touch vertical velocity < 0
    Jacobian_td = J_SwingConstraint_m(xT, p_full_T);
    vT_swingFoot = Jacobian_td * dqT;
    ineq = vT_swingFoot(2);

    % Bounds
    lower_bound = -inf;
    upper_bound = 0;
    name_bound = 'touch_down_vel_y';

    % Extend
    h = [h; ineq];
    lb = [lb; lower_bound];
    ub = [ub; upper_bound];
    boundName{end+1,1} = name_bound;

    %% 3. Lambda constraints
    % Bounds
    lower_bound = 1e-6;
    upper_bound = inf;
    name_bound = 'lambda_y';

    % Extend
    for k = 1:n_columns
        q_k = x(1:nq, k);
        dq_k = x(nq+1:end, k);
        u_k = u(:, k);
        p_full_k = p_full(:, k);
        lambda = compute_lambda([q_k;dq_k], u_k, p_full_k);
        % lambda = compute_lambda(q_k,dq_k, u_k, p_full_k);
        h = [h; lambda(2)];
        lb = [lb; lower_bound]; 
        ub = [ub; upper_bound];
        boundName{end+1,1} = name_bound;
    end                      
end
