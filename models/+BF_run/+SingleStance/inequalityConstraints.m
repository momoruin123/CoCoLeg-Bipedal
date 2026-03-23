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
    import BF_run.SingleStance.*;  % Import stance phase functions 

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
    
    %% Ground Clearance Constraint   
    % Swing foot pos. >= 0
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
    local_lower_bound = 0.8*0.06;

    for k = 1:(n_columns-1)
        lower_bound(k) = local_lower_bound;
    end  

    % Extend
    h = [h; ineq];
    lb = [lb; lower_bound];
    ub = [ub; upper_bound];
    for i = 1:n_columns
        boundName{end+1,1} = name_bound;
    end

    %% Lambda constraints
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
        h = [h; lambda(2)];
        lb = [lb; lower_bound]; 
        ub = [ub; upper_bound];
        boundName{end+1,1} = name_bound;
    end                      
end
