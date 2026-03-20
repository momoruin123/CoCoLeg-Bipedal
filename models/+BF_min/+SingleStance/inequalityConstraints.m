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

    %% Auxiliary Variable Bounds
    % Time step bounds
    a_lb.dt = 0.002;    % Minimum time step [s]
    a_ub.dt = 0.05;

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

    % Bounds
    lower_bound = 0;
    upper_bound = inf;
    name_bound = 'ground_clearance';

    % Extend
    h = [h; ineq];
    lb = [lb; lower_bound * ones(n_columns, 1)];
    ub = [ub; upper_bound * ones(n_columns, 1)];
    for i = 1:n_columns
        boundName{end+1,1} = name_bound;
    end

    %% Lift foot up limit
    % mid_k = floor(n_columns/2) + 1;
    % ineq = swingFoot_pos_y(mid_k);
    % 
    % % Bounds
    % lower_bound = 0.1;
    % upper_bound = inf;
    % name_bound = 'mid_step_lift_height';
    % 
    % % Extend
    % h = [h; ineq];
    % lb = [lb; lower_bound];
    % ub = [ub; upper_bound];
    % boundName{end+1,1} = name_bound;

    max_lift_height = 0.1;
    offset = floor(n_columns/4);
    for k = offset:n_columns-offset
        y_k = swingFoot_pos_y(k);
        t_swing = k/n_columns;
        target_y = max_lift_height * sin(pi * t_swing);

        ineq = y_k;
        lower_bound = target_y;
        upper_bound = inf;
        name_bound = sprintf('y_swing_phase_%d', k);

        h = [h; ineq];
        lb = [lb; lower_bound];
        ub = [ub; upper_bound];
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

%% Calculate the contact force between the ground and the stance foot
% function lambda = compute_lambda(q, dq, u, p)
%     import BF_min.SingleStance.*
%     temp = [0; 0];
%     q = vertcat(temp, q);
%     dq = vertcat(temp, dq);
%     nq = numel(q);
%     X_f = [q;dq];
% 
%     M_f = massMatrix_f(X_f, p);  
%     F_f= forces_f(X_f, u, p);
% 
%     J_Stance = J_StanceConstraint_f(X_f, p);
%     dJ_Stance = dJ_StanceConstraint_f(X_f, p);
% 
%     % gamma = dJ* dq
%     gamma = dJ_Stance * dq;
% 
%     A = [M_f    , -J_Stance';
%         J_Stance, zeros(size(J_Stance,1))];
%     B = [F_f;
%         -gamma];
% 
%     sol = A \ B;
%     lambda = sol(nq+1:end);
% end

%% function [h, lb, ub, boundName] = inequalityConstraints(config, Z, P_sym, N)
% % INEQUALITYCONSTRAINTS 定义非线性不等式约束（不含变量边界）
% %
% %   注意：变量的上下界（如 q, dq, u 的限制）请在外部通过 lbx/ubx 传入，
% %   不要写在这里，否则会极大降低 IPOPT 的收敛速度。
% 
%     import casadi.*
%     import BF_min.SingleStance.*;
% 
%     % --- 1. 基础尺寸获取 ---
%     [nx, nq, nu, ~] = get_sizes();
%     % na 包含优化参数和时间步长 dt
%     na = numel(config.optParameterNames) + config.optimizeTimeFlag;
%     n_rows = nx + nu + na;
%     n_columns = N+1; % 假设 Z 对应 N 个点
% 
%     % --- 2. 变量提取 ---
%     Z_mat = reshape(Z, n_rows, n_columns);
%     x = Z_mat(1:nx, :);
%     u = Z_mat(nx+1:nx+nu, :);
%     % q = x(1:nq, :);
%     % dq = x(nq+1:end, :);
% 
%     % --- 3. 初始化约束容器 ---
%     h = [];
%     lb = [];
%     ub = [];
%     boundName = {};
% 
%     %% A. 摆动腿地面清理 (Ground Clearance)
%     % 约束所有中间时刻摆动腿高度 >= 0
%     swingFoot_pos = cons_SwingFoot_m(x, P_sym); % 假设该函数支持矩阵输入
%     swing_y = swingFoot_pos(2, :);
% 
%     h = [h; swing_y'];
%     lb = [lb; zeros(n_columns, 1)];
%     ub = [ub; inf(n_columns, 1)];
%     for i = 1:n_columns
%         boundName{end+1} = sprintf('ground_clearance', i);
%     end
% 
%     %% B. 跨步中点抬腿高度 (Mid-step Lift Height)
%     % 强制要求步态中间点抬高至少 0.1m，避免拖地
%     mid_k = floor(n_columns/2) + 1;
%     swing_y_mid = swing_y(mid_k);
% 
%     h = [h; swing_y_mid];
%     lb = [lb; 0.1]; 
%     ub = [ub; inf];
%     boundName{end+1} = 'mid_step_lift_height';
% 
%     %% C. 落地与起跳速度方向 (Swing Foot Velocity)
%     % 起跳瞬间 (k=1) 垂直速度 > 0
%     J0 = J_SwingConstraint_m(x(:,1), P_sym);
%     v0 = J0 * x(nq+1:end, 1);
%     h = [h; v0(2)];
%     lb = [lb; 0];
%     ub = [ub; inf];
%     boundName{end+1} = 'lift_off_vel_y';
% 
%     % 落地瞬间 (k=N) 垂直速度 < 0
%     JT = J_SwingConstraint_m(x(:,end), P_sym);
%     vT = JT * x(nq+1:end, end);
%     h = [h; vT(2)];
%     lb = [lb; -inf];
%     ub = [ub; 0];
%     boundName{end+1} = 'touch_down_vel_y';
% 
%     %% D. 支撑腿不离地约束 (Normal Contact Force Lambda)
%     % 支撑腿垂直方向压力 lambda_y 必须为正 (防止跳起来)
%     % 注意：这步计算开销较大，如果收敛慢可以先注释掉
%     for k = [1, mid_k, n_columns] % 可以选几个关键点，不一定要全选
%         q_k = x(1:nq, k);
%         dq_k = x(nq+1:end, k);
%         u_k = u(:, k);
% 
%         lambda = compute_lambda(q_k, dq_k, u_k, P_sym);
%         h = [h; lambda(2)];
%         lb = [lb; 1e-6]; % 给予一个微小的正值
%         ub = [ub; inf];
%         boundName{end+1} = sprintf('lambda_y', k);
%     end
% 
%     % --- 4. 转换为 CasADi 表达式 ---
%     h = vertcat(h{:}); 
% 
% end
% 
% %% Calculate the contact force between the ground and the stance foot
% function lambda = compute_lambda(q, dq, u, p)
%     import BF_min.SingleStance.*
%     temp = [0; 0];
%     q = vertcat(temp, q);
%     dq = vertcat(temp, dq);
%     nq = numel(q);
%     X = [q;dq];
% 
%     M_f = massMatrix_f(X, p);  
%     F_f= forces_f(X, u, p);
% 
%     J_Stance = J_StanceConstraint_f(X, p);
%     dJ_Stance = dJ_StanceConstraint_f(X, p);
% 
%     % gamma = dJ* dq
%     gamma = dJ_Stance * dq;
% 
%     A = [M_f    , -J_Stance';
%         J_Stance, zeros(size(J_Stance,1))];
%     B = [F_f;
%         -gamma];
% 
%     sol = A \ B;
%     lambda = sol(nq+1:end);
% end 

%% function [h, lb, ub, boundName] = inequalityConstraints(config, Z, P_sym, N)
% % INEQUALITYCONSTRAINTS Define inequality constraints for stance phase optimization
% %
% %   Inputs:
% %     config - system configuration structure
% %     Z      - decision variables vector
% %     P_sym  - symbolic model parameters
% %     N      - number of collocation points
% %
% %   Outputs:
% %     h        - inequality constraint expressions in CASADI Symblic format
% %     lb, ub   - lower and upper bounds for constraints
% %     boundName- names of each constraint element
% %
% %   Author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025
% 
%     import casadi.*
%     import BF_min.SingleStance.*;  % Import stance phase functions
% 
%     % Load event function for contact force calculation
%     [eventFcn, ~, ~] = S2S.createEventJump();     
% 
%     % Determine number of auxiliary variables
%     na = numel(config.optParameterNames) + config.optimizeTimeFlag;
% 
%     % Get model dimensions
%     [nx, nq, nu, ~] = get_sizes();
% 
%     % Initialize additional constraint counters
%     n_add = 0;
%     add_lb = [];
%     add_ub = [];
%     add_names = {};
% 
%     % Reshape decision variables into matrix form
%     n_rows = nx + nu + na;
%     n_columns = numel(Z) / n_rows;
%     assert(mod(n_columns, 1) == 0 && n_columns > 0, 'n_columns must be a positive integer.');
% 
%     % Load numerical parameter values
%     P_num = getModelParameters(config, 0, 1);
% 
%     % Initialize parameter bound structures
%     paramValues = struct();
%     param_lb = struct();
%     param_ub = struct();
% 
%     % Get all parameter names from model
%     getParameterNames = str2func([config.model_name, '.getParameterNames']);
%     allParamNames = getParameterNames(config);
% 
%     % Set parameter values and bounds
%     for i = 1:numel(allParamNames)
%         paramName = allParamNames{i};
%         isOptimized = ismember(paramName, config.optParameterNames);
% 
%         % Check if parameter has inequality bounds override
%         hasIneqBound = isfield(config, 'inequalityParamBoundName') && ...
%                        ismember(paramName, config.inequalityParamBoundName);
% 
%         if isOptimized
%             % Use optimization bounds for optimized parameters
%             idxOpt = find(strcmp(paramName, config.optParameterNames));
% 
%             if hasIneqBound
%                 % Override with inequality bounds
%                 idxIneq = find(strcmp(paramName, config.inequalityParamBoundName));
%                 param_lb.(paramName) = config.inequalityParamLowerBound(idxIneq);
%                 param_ub.(paramName) = config.inequalityParamUpperBound(idxIneq);
%             else
%                 % Use standard optimization bounds
%                 param_lb.(paramName) = config.optParameterLowerBound(idxOpt);
%                 param_ub.(paramName) = config.optParameterUpperBound(idxOpt);
%             end
% 
%             paramValues.(paramName) = config.optParameterLowerBound(idxOpt);
%         else
%             % Use fixed values for non-optimized parameters
%             idx = getParameterIndex(paramName, config);
%             paramValues.(paramName) = P_num(idx);
% 
%             if hasIneqBound
%                 % Apply inequality bounds even to fixed parameters
%                 idxIneq = find(strcmp(paramName, config.inequalityParamBoundName));
%                 param_lb.(paramName) = config.inequalityParamLowerBound(idxIneq);
%                 param_ub.(paramName) = config.inequalityParamUpperBound(idxIneq);
%             else
%                 % Fixed parameters have equal lower/upper bounds
%                 param_lb.(paramName) = P_num(idx);
%                 param_ub.(paramName) = P_num(idx);
%             end
%         end
%     end
% 
%     % % Calculate natural frequency bounds for time step constraints
%     % eigenfreq_min = 1/(2*pi) * sqrt(param_lb.k_l / paramValues.m);  % [Hz]
%     % eigenfreq_max = 1/(2*pi) * sqrt(param_ub.k_l / paramValues.m);  % [Hz]
%     % T_max = 1 / eigenfreq_min; % Maximum period [s]
%     % T_min = 1 / eigenfreq_max; % Minimum period [s]
% 
%     % Initialize constraint arrays
%     h = MX.zeros(2 * numel(Z));
%     idx_h = 0;
% 
%     %% State Bounds
%     % State bounds
%     state_lb = config.state_lb;
%     state_ub = config.state_ub;
% 
%     x_lb.phi = state_lb(1);
%     x_ub.phi = state_ub(1);
% 
%     x_lb.alphaL = state_lb(2);
%     x_ub.alphaL = state_ub(2);
% 
%     x_lb.alphaR = state_lb(3);
%     x_ub.alphaR = state_ub(3);
% 
%     x_lb.betaL = state_lb(4);
%     x_ub.betaL = state_ub(4);
% 
%     x_lb.betaR = state_lb(5);
%     x_ub.betaR = state_ub(5);
% 
%     %% Input Bounds
%     input_lb = config.input_lb;
%     input_ub = config.input_ub;
% 
%     u_lb.u_hl = input_lb(1);  
%     u_ub.u_hl = input_ub(1);  
% 
%     u_lb.u_hr = input_lb(2);  
%     u_ub.u_hr = input_ub(2);
% 
%     u_lb.u_kl = input_lb(3);  
%     u_ub.u_kl = input_ub(3);
% 
%     u_lb.u_kr = input_lb(4);  
%     u_ub.u_kr = input_ub(4);
% 
%     %% Auxiliary Variable Bounds
%     % Time step bounds
%     epsilon = 1e-2 / N;  % Minimum time step [s]
%     a_lb.dt = epsilon;
%     % =====================================================
%     % ==================== NEED ADJUST ====================
%     % =====================================================
%     % a_ub.dt = 1 / N;  % Maximum time step
%     a_ub.dt = 0.01;
% 
%     % Optimization parameter bounds
%     for i = 1:numel(config.optParameterNames)
%         paramName = config.optParameterNames{i};
%         a_lb.(paramName) = config.optParameterLowerBound(i);
%         a_ub.(paramName) = config.optParameterUpperBound(i);
%     end
% 
%     %% Additional Dynamic Constraints
%     % Extract variables from decision vector
%     Z = reshape(Z, n_rows, n_columns);
%     Z_extended = Z;
%     x = Z(1:nx, :);
%     u = Z(nx+1:nx+nu, :);
%     a = Z(nx+nu+1:n_rows, :);
% 
%     %% Ground Clearance Constraint   
%     % Swing foot >=0
%     % epsilon_foot = 0;  % Foot clearance tolerance
%     % swingFoot_pos = cons_SwingFoot_m(x, P_sym);
%     % swingFoot_pos_y = swingFoot_pos(2,:);
%     % 
%     % add_lb.ground_clearance = -epsilon_foot;  % Minimum foot height
%     % add_ub.ground_clearance = inf;   
%     % 
%     % Z_extended = [Z_extended; swingFoot_pos_y];
%     % n_add = n_add + 1;
% 
%     %% Lambda constraints
%     % epsilon_lambda = 1e-6;
%     % q = x(1:nq,:);
%     % dq = x(nq+1:end,:);
%     % lambda = compute_lambda(q, dq, u, P_sym);
%     % lambda_y = lambda(2);
%     % 
%     % % contact force lambda of stance foot must be positive
%     % add_lb.contact_lambda = lambda_y - epsilon_lambda;
%     % add_ub.contact_lambda = inf;                          
%     % 
%     % Z_extended = [Z_extended; lambda_y];
%     % n_add = n_add + 1;
% 
%     %% 1. Lift foot up limit
%     epsilon_liftup = 0.1;
%     mid_k = floor(N/2) + 1;
%     x_k = x(:,mid_k);
%     swingFoot_pos_m = cons_SwingFoot_m(x_k, P_sym);
%     swingFoot_pos_y_m = swingFoot_pos_m(2);
% 
%     % Minimum foot height in the middle
%     % add_lb.liftup_midpos = epsilon_liftup;  
%     % add_ub.liftup_midpos = inf;   
% 
%     % Z_extended = [Z_extended; swingFoot_pos_y_m];
%     % n_add = n_add + 1;
%     %% 2. 3. Swing foot velocity constraints
%     % extract initial state and final state
%     x0 = x(:,1);                % Initial state
%     q0 = x0(1:nq);              
%     dq0 = x0(nq+1:end);         
%     xT = x(:,end);              % final state
%     qT = xT(1:nq);       
%     dqT = xT(nq+1:end);  
% 
%     % Lift off: move upwards
%     Jacobian_lift = J_SwingConstraint_m(x0, P_sym);
%     vel_swingFoot_lift = Jacobian_lift * dq0;
%     lift_vel_y = vel_swingFoot_lift(2);      % Lift vertical velocity > 0
% 
%     % Touch down: move downwards
%     Jacobian_touch = J_SwingConstraint_m(xT, P_sym);
%     vel_swingFoot_touch = Jacobian_touch * dqT;
%     touch_vel_y = vel_swingFoot_touch(2);   % Touch vertical velocity < 0
% 
%     % Extend variable set with lift and touch velocity
%     % Z_extended = [Z_extended; lift_vel]; 
%     % Z_extended = [Z_extended; touch_vel];        
%     % n_add = n_add + 2;  
% 
%     %% putting it all together
%     n_Z_extended = n_columns * (n_rows + n_add);
%     Z_extended = reshape(Z_extended, n_Z_extended, 1);
% 
%     % Convert bound structures to vectors
%     [x_lb, x_lb_names] = s2v(x_lb);
%     [x_ub, ~] = s2v(x_ub);
%     [u_lb, u_lb_names] = s2v(u_lb);
%     [u_ub, ~] = s2v(u_ub);
%     [a_lb, a_lb_names] = s2v(a_lb);
%     [a_ub, ~] = s2v(a_ub);
%     % [add_lb, add_lb_names] = s2v(add_lb);
%     % [add_ub, ~] = s2v(add_ub);
% 
%     % Create bound vectors for all collocation points
%     % lb = repmat([x_lb; u_lb; a_lb(1:na); add_lb], n_columns, 1);
%     % ub = repmat([x_ub; u_ub; a_ub(1:na); add_ub], n_columns, 1);
%     % boundName = repmat([x_lb_names'; u_lb_names'; a_lb_names(1:na)'; add_lb_names'], n_columns, 1);
%     lb = repmat([x_lb; u_lb; a_lb(1:na)], n_columns, 1);
%     ub = repmat([x_ub; u_ub; a_ub(1:na)], n_columns, 1);
%     boundName = repmat([x_lb_names'; u_lb_names'; a_lb_names(1:na)'], n_columns, 1);
% 
% 
%     %%
%     % Remove unbounded constraints for efficiency
%     unnecessaryConst = (lb == -inf) & (ub == inf);
%     casadiIdx = find(unnecessaryConst == 0);
%     lb = lb(~unnecessaryConst);
%     ub = ub(~unnecessaryConst);
%     boundName = boundName(~unnecessaryConst);
% 
%     % Assign constraint expressions
%     h(idx_h+1:numel(lb)) = Z_extended(casadiIdx);
%     idx_h = idx_h + numel(lb);
% 
%     %% Add individual constraints
%     % % 1. Lift foot up minimal height > 0.1
%     % idx_h = idx_h + 1;
%     % h(idx_h) = swingFoot_pos_y_m;
%     % lb = vertcat(lb, 0.1);          
%     % ub = vertcat(ub, inf); 
%     % boundName = [boundName; {'min_lift_up_height'}]; 
%     % 
%     % % 2. Lift off velocity limit：lift_vel > 0
%     % idx_h = idx_h + 1;
%     % h(idx_h) = lift_vel_y;
%     % lb = vertcat(lb, 0);          
%     % ub = vertcat(ub, inf);        
%     % boundName = [boundName; {'swing_vel_lift_off'}]; 
%     % 
%     % % 3. Touch down velocity limit：-touch_vel > 0
%     % idx_h = idx_h + 1;
%     % h(idx_h) = -touch_vel_y;
%     % lb = vertcat(lb, 0);          
%     % ub = vertcat(ub, inf);        
%     % boundName = [boundName; {'swing_vel_touch_down'}];
% 
%     %% Truncate to actual number of constraints
%     h = h(1:idx_h);
%     boundName = boundName(1:idx_h);
% 
% end
% 
% %% Calculate the contact force between the ground and the stance foot
% function lambda = compute_lambda(q, dq, u, p)
%     import BF_min.SingleStance.*
%     temp = [0; 0];
%     q = vertcat(temp, q);
%     dq = vertcat(temp, dq);
%     nq = numel(q);
%     X = [q;dq];
% 
%     M_f = massMatrix_f(X, p);  
%     F_f= forces_f(X, u, p);
% 
%     J_Stance = J_StanceConstraint_f(X, p);
%     dJ_Stance = dJ_StanceConstraint_f(X, p);
% 
%     % gamma = dJ* dq
%     gamma = dJ_Stance * dq;
% 
%     A = [M_f    , -J_Stance';
%         J_Stance, zeros(size(J_Stance,1))];
%     B = [F_f;
%         -gamma];
% 
%     sol = A \ B;
%     lambda = sol(nq+1:end);
% end