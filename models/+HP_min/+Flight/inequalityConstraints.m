function [h, lb, ub, boundName] = inequalityConstraints(config, Z, P_sym, N)
    import casadi.*
    % Import all public members from mypackage
    import HP_min.Flight.*;

    na  = numel(config.optParameterNames) + config.optimizeTimeFlag;
    
    % find a way to get nx, nq, nu and np.
    [nx, nq, nu, np] = get_sizes();

    % Initialize additional constraint counters
    n_add = 0;
    add_lb = [];
    add_ub = [];
    
    % Reshape decision variables into matrix form
    n_rows = nx + nu + na;
    n_columns = numel(Z) / n_rows;
    assert(mod(n_columns, 1) == 0 && n_columns > 0, 'n_columns must be a positive integer.');

    % Load numerical parameter vector (e.g., from getFullParameters)
    P_num = getModelParameters(config, 0, 1);
    
    % Initialize structure to hold all parameter values
    paramValues = struct();

    % Loop over all known parameter names from the model
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

    % Initialize residuals array
    h     = MX.zeros(2*numel(Z));
    idx_h = 0;

    %% State bounds for flight phase
    % y-position
    x_lb.y     =  paramValues.r_f;
    x_ub.y     =  5*paramValues.l_0;
        
    % length (l)
    x_lb.l     =  1/2 * paramValues.l_0;
    x_ub.l     =  3/2 * paramValues.l_0;
    
    % y-velocity
    x_lb.dy    = -sqrt(2*paramValues.g*x_ub.y);
    x_ub.dy    = +sqrt(2*paramValues.g*x_ub.y);
    
    % Leg extension velocity bounds 
    x_lb.dl = -sqrt(5 * paramValues.g * x_ub.l);   
    x_ub.dl = sqrt(5 * paramValues.g * x_ub.l);

    %% input bounds in Flight
    epsilon_l     = 4;

    u_lb.u_l    = -epsilon_l*paramValues.m*paramValues.g; 
    u_ub.u_l    = +epsilon_l*paramValues.m*paramValues.g;

    %% auxilliary bounds
    epsilon   = 3e-4/N; %[s]
    a_lb.dt   = epsilon;
    a_ub.dt   = +4*T_max/N;

    % Loop over all optimization parameters
    for i = 1:numel(config.optParameterNames)
        paramName = config.optParameterNames{i};
        a_lb.(paramName) = config.optParameterLowerBound(i);
        a_ub.(paramName) = config.optParameterUpperBound(i);
    end

    %% Ground clearance in flight
    epsilon_foot = 3e-3;
    Z = reshape(Z, nx+nu+na, n_columns);
    x = Z(1:nx, :);
    foot_pos = x(1,:)-x(2,:);  % y-l

    add_lb.foot    = paramValues.r_f - epsilon_foot;
    add_ub.foot    = x_ub.y;

    Z_extended = [Z; foot_pos];
    Z_extended = reshape(Z_extended, n_columns*(nx+nu+na+1), 1);

    %% putting it all together

    [x_lb, x_lb_names] = s2v(x_lb);
    [x_ub, ~] = s2v(x_ub);
    [u_lb, u_lb_names] = s2v(u_lb);
    [u_ub, ~] = s2v(u_ub);
    [a_lb, a_lb_names] = s2v(a_lb);
    [a_ub, ~] = s2v(a_ub);
    [add_lb, add_lb_names] = s2v(add_lb);
    [add_ub, ~] = s2v(add_ub);

    lb = repmat([x_lb; u_lb; a_lb(1:na); add_lb], n_columns, 1);
    ub = repmat([x_ub; u_ub; a_ub(1:na); add_ub], n_columns, 1);

    boundName = repmat([x_lb_names'; u_lb_names'; a_lb_names(1:na)'; add_lb_names], n_columns, 1);

    % remove unnecessary constraints
    unnecessaryConst = (lb == -inf) & (ub == inf);
    casadiIdx = find(unnecessaryConst == 0);
    lb = lb(~unnecessaryConst);
    ub = ub(~unnecessaryConst);
    boundName = boundName(~unnecessaryConst);

    h(idx_h+1:numel(lb)) = Z_extended(casadiIdx);
    idx_h = idx_h + numel(lb);

    %%
    % truncate h 
    h = h(1:idx_h);
    boundName = boundName(1:idx_h);

end
