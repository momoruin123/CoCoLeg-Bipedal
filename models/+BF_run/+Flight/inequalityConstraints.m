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
    import BF_run.Flight.*;  % Import stance phase functions 

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
    swingFoot_pos = pos_SwingFoot_f(x, p_full);
    swingFoot_pos_y = swingFoot_pos(2,:)';
    ineq = swingFoot_pos_y;

    % Bounds
    lower_bound = 0;
    upper_bound = inf;
    name_bound = 'swing_foot_ground_clearance';

    % Extend
    h = [h; ineq];
    lb = [lb; lower_bound * ones(n_columns, 1)];
    ub = [ub; upper_bound * ones(n_columns, 1)];
    for i = 1:n_columns
        boundName{end+1,1} = name_bound;
    end

    % Stance foot pos. >= 0
    stanceFoot_pos_y = x(2,:)';
    ineq = stanceFoot_pos_y;

    % Bounds
    lower_bound = 0;
    upper_bound = inf;
    name_bound = 'stance_foot_ground_clearance';

    % Extend
    h = [h; ineq];
    lb = [lb; lower_bound * ones(n_columns, 1)];
    ub = [ub; upper_bound * ones(n_columns, 1)];
    for i = 1:n_columns
        boundName{end+1,1} = name_bound;
    end

    %% Vertical clearance constraint betweent two foot
    % At the last segment of flight phase, stance foot should be higher
    % than swing foot (from last SingleStance phase)

    % Extract necessary state
    xT  = x(:, end);
    p_full_T = p_full(:, end);

    % Vertical stance foot - Vertical swing foot >= 0 at the end segment
    stanceFoot_pos_y = xT(2,end);
    swingFoot_pos    = pos_SwingFoot_f(xT, p_full_T);
    swingFoot_pos_y  = swingFoot_pos(2,end)';
    ineq = stanceFoot_pos_y - swingFoot_pos_y;

    % Bounds
    lower_bound = 0;
    upper_bound = inf;
    name_bound = 'bipedal_foot_clearance';

    % Extend
    h  = [h ; ineq];
    lb = [lb; lower_bound];
    ub = [ub; upper_bound];
    boundName{end+1,1} = name_bound;

    %% Swing foot velocity constraints
    % Extract necessary state
    dqT  = dq(:, end);

    % Touch down: Touch vertical velocity < 0
    Jacobian_td  = J_SwingConstraint_f(xT, p_full_T);
    vT_swingFoot = Jacobian_td * dqT;
    ineq = -vT_swingFoot(2);

    % Bounds
    lower_bound = 0;
    upper_bound = inf;
    name_bound = 'touch_down_vel_y';

    % Extend
    h  = [h; ineq];
    lb = [lb; lower_bound];
    ub = [ub; upper_bound];
    boundName{end+1,1} = name_bound;

end
