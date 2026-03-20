function [h, lb, ub, boundName] = inequalityConstraints(config, Z, P_sym, N)
    import casadi.*
    % Import all public members from mypackage
    import HP_min.Stance.*;

    na  = numel(config.optParameterNames) + config.optimizeTimeFlag;
    % find a way to get nx, nq, nu and np.
    [nx, nq, nu, np] = get_sizes();

    % Initialize additional constraint counters
    n_add = 0;
    add_lb = [];
    add_ub = [];

    n_rows = nx+nu+na;
    n_columns = numel(Z)/n_rows;
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
    h = MX.zeros(2*numel(Z));
    idx_h = 0;

    %% State bounds in stance phase

    % length (l)
    x_lb.l     =  1/2 * paramValues.l_0;
    x_ub.l     =  3/2 * paramValues.l_0;
   
    % Leg extension velocity bounds 
    x_lb.dl = -sqrt(5 * paramValues.g * x_ub.l);   
    x_ub.dl = sqrt(5 * paramValues.g * x_ub.l);

    % Alternative to be used in a FSSF config to limit the stance to a single swing
    % epsilon_dl = (1/12)*sqrt(2*paramValues.g*x_ub.l);
    % x_lb.dl    = config.ineqConstraints.coef(1)*(sqrt(2*paramValues.g*x_ub.l))-epsilon_dl;   
    % x_ub.dl    = config.ineqConstraints.coef(2)*(sqrt(2*paramValues.g*x_ub.l))+epsilon_dl;   

    %% input bounds in Stance
    epsilon_l = 4;
    u_lb.u_l  = -epsilon_l*paramValues.m*paramValues.g;  
    u_ub.u_l  = +epsilon_l*paramValues.m*paramValues.g;
    
    %% auxilliary bounds
    epsilon   = 3e-4/N; %[s]
    a_lb.dt_Stance   = epsilon;
    a_ub.dt_Stance   = +4*T_max/N;

    % Loop over all optimization parameters
    for i = 1:numel(config.optParameterNames)
        paramName = config.optParameterNames{i};
        a_lb.(paramName) = config.optParameterLowerBound(i);
        a_ub.(paramName) = config.optParameterUpperBound(i);
    end

    %% Positive contact force in stance
    % Extract stance phase variables
    Z = reshape(Z, n_rows, n_columns);
    x = Z(1:nx, :);
    u = Z(nx+1:nx+nu, :);
    a = Z(nx+nu+1:nx+nu+na, :);

    % we know we can find the contact force in the event Stance to flight 
    [eventFcn, ~, ~] = HP_min.Stance.S2F.createEventJump();     
    
    p_matrix   = repmat(P_sym, 1, size(a,2));
    p_aux      = a(end-na+2:end, :);
    p_full     = [p_matrix; p_aux];  
    pos_lambda = eventFcn(x, u, p_full);

    add_lb.lambda    = -1e-8;
    add_ub.lambda    = +inf;

    Z_extended = [Z; pos_lambda];
    n_add = n_add + 1;

    %% One swing in Stance
    % Alternative way to be used in a FSSF config to limit the stance to a single swing
    epsilon_half_swing = 7/12;
    a_dt = a(end-na+1, :);
    % eigenfreq = 1/(2*pi)*sqrt(p_full(getParameterIndex('k_l', config), :)/paramValues.m);
    half_swing = a_dt - epsilon_half_swing*(1/eigenfreq_max)/N;

    add_lb.half_swing = -inf;
    add_ub.half_swing = 0;

    Z_extended = [Z_extended; half_swing];
    n_add = n_add + 1;

    %% limit leg displacement to natural 
    l_multiple_k = (x(1,:)-paramValues.l_0).*p_full(getParameterIndex('k_l', config), :);
    add_lb.l_multiple_k = -paramValues.m*paramValues.g -u_ub.u_l;  %u_ub.u_l-2*paramValues.m*paramValues.g;
    add_ub.l_multiple_k = -paramValues.m*paramValues.g +u_ub.u_l;      %u_ub.u_l;

    Z_extended = [Z_extended; l_multiple_k];
    n_add = n_add + 1;
    %% putting it all together
    n_Z_extended = n_columns*(n_rows + n_add);
    Z_extended = reshape(Z_extended, n_Z_extended, 1);

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

    boundName = repmat([x_lb_names'; u_lb_names'; a_lb_names(1:na)'; add_lb_names'], n_columns, 1);

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