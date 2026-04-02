%% Finding Good Starting Points for Hopping Robot Optimization


% This file helps find good starting guesses for the hopping robot trajectory optimization. 
% It serves solely as an example of how feasible seed points can be generated. 
% The presented approach is computationally intensive, as it explores multiple trials 
% and retains those that yield valid results rather than relying on a well-defined heuristic.

% Here's what we do step by step:

% We start by creating a big list of different starting conditions:
% - Various leg angles at launch
% - Different hopping speeds (forward and upward)
% - Different leg spring stiffnesses
% - Different control force levels

% For each combination in our list:
% 1. We simulate the robot hopping forward (using impulsive forces, PD controller)
% 2. We check if the hop looks physically reasonable
% 3. We measure how fast the robot moves on average
% 4. We find which target speed (from our list) is closest
% 5. We run an optimization without optimizing the leg stiffness, at the average speed estiamted before
% 6. We run a second optimization to adjust leg stiffness using the solution of the first as an initial guess
% 7. We check if optimization worked well and save it if it's the best at
% that grid point (k_l, v_avg)

% At the end, we have a collection of good starting points. We use these
% starting points as seeds in the expolaration (see main_optimization)

% Warning: This approach is tailored specifically to the hopping model in minimal coordinates and would need significant changes for other types of robots.

% author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

%%
% Generate dated filename for results
dated_filename = append_date('WS');

% Load model configuration
configfile_BF_min;

% Save original configuration and set backup parameters
config_0 = config;
backup_batchSize = 50;
% backup_batchSize = 1000;

% Set hip compliance parameter
config.paramValues.k_h = 20/256;
config.paramValues.k_k = 20/256;

%% Initialize Model Parameters and State
% Get model parameters and initial state
p_vec = getModelParameters(config, 0, 0);
X0    = config.x0;

%% Initialize Simulation Settings
ode.settings = odeset('RelTol', 1e-7, 'AbsTol', 1e-8);
ode.solver   = @ode45;
dt = 3e-2;  % Time step

% Define parameters to optimize
optParameterNames = config.optParameterNames;

%% Define Parameter Grids for Warm Start
% grid sample of velocity
v_avg_grid   = [0.4, 0.7];
% v_avg_grid   = [0.01, 0.1:0.1:1.5];
% v_avg_grid   = 0.001:0.1:1.501;
num_v_avg    = numel(v_avg_grid);

% Latin hypercube sample of params
% k_h_grid, k_k_grid, L_grid, H_grid, phi_grid

% Config. of sampling
n_samples = 1;
lb = [0, 0, 0.4, 0.4, -pi/3]; % lower bounds
ub = [100, 100, 1, 0.8, pi/3]; % upper bounds

% lb = [0/256, 0/256, 0.6, 0.5, 0]; % lower bounds
% ub = [0/256, 0/256, 0.6, 0.5, 0]; % upper bounds
n_g = numel(lb);
n_p = numel(optParameterNames);    % number of optimized params

% Latin_hypercube_sampling
input_matrix = latin_hypercube_sampling(lb, ub, n_samples);
k_h_grid   = input_matrix(:,1);     % hip stiffness
k_k_grid   = input_matrix(:,2);     % knee stiffness
L_grid     = input_matrix(:,3);     % step length
H_grid     = input_matrix(:,4);     % hip height
phi_grid   = input_matrix(:,5);     % pitch angle

% Flatten grids into vectors
input_matrix = repmat(input_matrix, num_v_avg, 1);

% Define parameter grids for warm start exploration
k_h_vector   = input_matrix(:,1);     % hip stiffness
k_k_vector   = input_matrix(:,2);     % knee stiffness
L_vector     = input_matrix(:,3);     % step length
H_vector     = input_matrix(:,4);     % hip height
phi_vector   = input_matrix(:,5);     % pitch angle
v_avg_vector = repelem(v_avg_grid', n_samples);     % velocity

%% Initialize Data Storage
% Test trajectory to determine state dimension
trajTest = flow(config, dt);
for j = 1:numel(trajTest)
    trajTest(j).p = [];
end
config.optParameterNames = {};
Z_test   = interpolate_traj2Z(config, trajTest);
n_z = numel(Z_test);

% Initialize result arrays
operatingPoint_array = v_avg_grid';             % Save all avg_vel point
cost_array = inf(num_v_avg, 1);                 % Save optimal cost for every avg_vel point
Z_array = zeros(num_v_avg, n_z);                % Optimal state vector Z
warmStart_array = nan(num_v_avg, n_g);          % Optimal warm-start parameters
gridPoint_array = nan(num_v_avg, n_p+1);        % 
numConverged_array = zeros(num_v_avg, 1);       % Number of convergences
largestNorm     = zeros(num_v_avg, 1);          % Maximum norm of constraint violation for every avg_vel point
numberTrials    = zeros(num_v_avg, 1);          % Warm start attempts for every avg_vel point

% Initialize tracking arrays
trackParams = nan(n_g+1, numel(v_avg_vector));          
trackCost = nan(2, numel(v_avg_vector));
trackOptParam = nan(2, n_p, numel(v_avg_vector)); % 

%% Restart the parallel processing pool. 
disp('starting warm start')
% delete(gcp('nocreate'));  % close any existing pool
% setenv('OMP_NUM_THREADS', '1');
% setenv('MKL_NUM_THREADS', '1');
% parpool('local');  % Number of local core

%% Main Optimization Loop
kk = 0;
failedIterations = [];

while kk < numel(v_avg_vector)
    disp(kk);
    % Preallocate results structure for current batch
    results = repmat(struct( ...
        'idx', [], ...
        'success', false, ...
        'cost', inf, ...
        'Z', [], ...
        'operatingPoint', nan, ...
        'gridPoint', nan(1,2), ...
        'warmStart', nan(1,4), ...
        'trackCost', nan(2,1), ...
        'trackOptParam', nan(2, n_p), ...
        'trackParams', nan, ...
        'normToOptimal', 0, ...
        'failReason', "" ... 
    ), 1, min(kk+backup_batchSize, numel(v_avg_vector))-kk);

    % Parallel processing of batch
    for i = 1:min(kk+backup_batchSize, numel(v_avg_vector))-kk
        localResult = results(i);

        % try
            % Configure current iteration
            config_i = config;
            config_i.paramValues.k_h = k_h_vector(i+kk);
            config_i.paramValues.k_k = k_k_vector(i+kk);
            config_i.optParameterNames = {};
            p_vec = getModelParameters(config_i, 0, 0);

            % Set initial conditions and control inputs
            config_i.x0 = Bipedal_IK(config, L_vector(i+kk), H_vector(i+kk), phi_vector(i+kk), v_avg_vector(i+kk));
            % config_i.x0 = X0;

            % Cheak physical feasibility
            if isempty(config_i.x0)
                localResult.failReason = 'IK infeasible';
                results(i) = localResult;
                continue;
            end

            % Generate initial guess trajectory
            % 1. Interpolation
            T_i_guess = L_vector(i+kk) / v_avg_vector(i+kk);
            trajInit = initial_guess(config_i, T_i_guess, 0);

            traj_full = trajToFulltraj(config, trajInit, 0);
            getAnimationRABBIT(traj_full, 0, config,[]);


            Z_init = interpolate_traj2Z(config_i, trajInit);

            % 2. Simulation
            % trajInit = flow(config, dt);

            % =========================NEED ADJUST========================
            % Find closest target velocity
            [~, idx] = min(abs(v_avg_grid - v_avg_vector(i+kk)));
            localResult.idx = idx;
            localResult.trackParams = [v_avg_grid(idx); k_h_vector(i+kk); k_k_vector(i+kk); L_vector(i+kk); H_vector(i+kk); phi_vector(i+kk)];
            localResult.operatingPoint = v_avg_grid(idx);

            % Step 1: Optimize with fixed stiffness
            config_i.operatingCond.v_avg  = v_avg_vector(i+kk);
            config_i.optConfig.tolerance = 1;
            [Z_fixed, P_fixed, w_fixed, r_fixed, ~, h_fixed, bounds_fixed, ~] = getConstraintsAndCosts(config_i, trajInit);

            [solutionIPOPT, stats] = solveCasadi(config_i, Z_fixed, P_fixed, r_fixed, h_fixed, w_fixed, bounds_fixed, trajInit);

            if ~strcmp(stats.return_status, 'Solve_Succeeded')
                localResult.failReason = 'opt 1 failed';
                results(i) = localResult;
                continue;
            end

            solution = full(solutionIPOPT.x);
            localResult.trackCost(1) = full(solutionIPOPT.f);
            localResult.trackOptParam(1, 1) = k_h_vector(i+kk);
            localResult.trackOptParam(1, 2) = k_k_vector(i+kk);

            % Interpolate solution
            trajSol = interpolate_Z2traj(config_i, solution);

            % Step 2: Optimize with adaptive stiffness
            % current p guess
            current_p_guess = [k_h_vector(i+kk), k_k_vector(i+kk)];
            p_matrix = repmat(current_p_guess, size(trajSol(1).x,1), 1);

            [trajSol.p] = deal(p_matrix);
            
            config_i.optParameterNames = optParameterNames;
            config_i.optParameterInit = [config_i.paramValues.k_h, config_i.paramValues.k_k];
            config_i.optConfig.tolerance = 2;
            [Z_adaptive, P_adaptive, w_adaptive, r_adaptive, ~, h_adaptive, bounds_adaptive, names] = getConstraintsAndCosts(config_i, trajSol);

            [solutionIPOPT, stats] = solveCasadi(config_i, Z_adaptive, P_adaptive, r_adaptive, h_adaptive, w_adaptive, bounds_adaptive, trajSol);

            if ~strcmp(stats.return_status, 'Solve_Succeeded')
                localResult.failReason = 'opt 2 failed';
                results(i) = localResult;
                continue;
            end

            solution = full(solutionIPOPT.x);
            traj     = interpolate_Z2traj(config_i, solution);
            localResult.trackCost(2) = full(solutionIPOPT.f);
            localResult.trackOptParam(2, :) = traj(1).p(1,:);

            % Check if fixed stiffness performed better
            if localResult.trackCost(1) < localResult.trackCost(2)
                localResult.failReason = 'fixed stiffness has better cost';
                results(i) = localResult;
                continue;
            end

            % Store successful result
            [traj.p] = deal([]);
            localResult.Z             = interpolate_traj2Z(config_i, traj);
            localResult.cost          = localResult.trackCost(2);
            localResult.gridPoint     = [v_avg_grid(idx), localResult.trackOptParam(2, :)];
            localResult.warmStart     = [k_h_vector(i+kk), k_k_vector(i+kk), L_vector(i+kk), H_vector(i+kk), phi_vector(i+kk)];
            localResult.normToOptimal = norm(Z_init-localResult.Z)^2;
            localResult.success       = true;

        % catch ME
        %     localResult.failReason = ['unhandled error: ', ME.identifier];
        %     results(i) = localResult;
        %     continue;
        % end
                    
        results(i) = localResult;

    end

    %%
    % Process results from current batch
    for i = 1:numel(results)
        % Update tracking arrays
        trackCost(:, i+kk)        = results(i).trackCost;
        trackOptParam(:, :, i+kk) = results(i).trackOptParam;
        trackParams(:, i+kk)      = results(i).trackParams;

        if ~results(i).success || numel(results(i).failReason)>1
            failedIterations = [failedIterations; {i+kk, results(i).failReason}];
            continue;
        elseif isempty(results(i).idx)
            failedIterations = [failedIterations; {i+kk, 'no idx found'}];
            continue;
        end
        
        idx = results(i).idx;
        numberTrials(idx) = numberTrials(idx) + 1;
        
        % Update best solution if current is better
        if results(i).success && results(i).cost < cost_array(idx)
            Z_array(idx,:)         = results(i).Z;
            cost_array(idx)        = results(i).cost;
            gridPoint_array(idx,:) = results(i).gridPoint;
            warmStart_array(idx,:) = results(i).warmStart;
            operatingPoint_array(idx,:) = results(i).operatingPoint;
            numConverged_array(idx) = 1;
            largestNorm(idx) = results(i).normToOptimal;
        % Handle similar solutions
        elseif norm(results(i).cost - cost_array(idx)) < 1e-2 && norm(Z_array(idx,:) - results(i).Z') < 1e-1 && norm(gridPoint_array(idx,2)-results(i).gridPoint(2)) < 1   
            numConverged_array(idx) = numConverged_array(idx)+ 1;
            if results(i).normToOptimal > largestNorm(idx)
                largestNorm(idx) = results(i).normToOptimal;
                warmStart_array(idx,:) = results(i).warmStart;
            end
        end

    end 
    
    % Save intermediate results
    save_data([dated_filename, '_UNFINISHED'], 'o', 'Z_array', Z_array, 'gridPoint_array', gridPoint_array, 'cost_array', cost_array, ...
        'warmStart_array', warmStart_array, 'numConverged_array', numConverged_array ,'largestNorm', largestNorm , 'numberTrials', numberTrials, 'config', config);
    
    % Move to next batch
    kk = kk+backup_batchSize;

end

%%

%% Save Final Results and Analysis
diagnosis = {};
% diagnosis.k_h_grid         = k_h_vector;
% diagnosis.k_k_grid         = k_k_vector;
% diagnosis.L_grid           = L_vector;
% diagnosis.H_grid           = H_vector;
% diagnosis.phi_grid         = phi_vector;
% diagnosis.v_avg_grid       = v_avg_vector;
diagnosis.k_h_grid         = k_h_grid;
diagnosis.k_k_grid         = k_k_grid;
diagnosis.L_grid           = L_grid;
diagnosis.H_grid           = H_grid;
diagnosis.phi_grid         = phi_grid;
diagnosis.v_avg_grid       = v_avg_grid;
diagnosis.failedIterations = failedIterations;
diagnosis.trackCost        = trackCost;
diagnosis.trackOptParam    = trackOptParam;
diagnosis.trackParams      = trackParams;
% diagnosis.trackParamMatrix = [input_matrix, v_avg_vector].';

%%
% Save final data
filename = save_data(dated_filename, 'v', 'Z_array', Z_array, 'gridPoint_array', gridPoint_array, 'cost_array', cost_array, ...
    'warmStart_array', warmStart_array, 'numConverged_array', numConverged_array ,'largestNorm', largestNorm , 'diagnosis', diagnosis, 'numberTrials', numberTrials, 'config', config_0);

% Analyze and document results
% analyzeFailuresWarmStart(failedIterations, numel(v_avg_vector), [filename,'.txt']);
% documentInequalityConstraints(config, [filename,'_ineq.txt'])

disp('warm start finished')

%% Latin_hypercube_sampling function
function input_samples = latin_hypercube_sampling(lb, ub, n_samples)
    n_params = numel(lb);

    % 0~1 sampling
    X_norm = lhsdesign(n_samples, n_params, 'criterion', 'maximin', 'iterations', 10);

    % Mapping from 0~1 to bounds
    input_samples = lb + X_norm .* (ub - lb);
    
end

%% 
% for i = 1:numel(v_avg_vector)
%     v_avg_val = v_avg_vector(i);
%     L_val     = L_vector(i);
%     H_val     = H_vector(i);
%     phi_val   = phi_vector(i);
%     k_h_val   = k_h_vector(i);
%     k_k_val   = k_k_vector(i);
% 
%     x0_i = Bipedal_IK(config, L_val, H_val, phi_val, v_avg_val);
% 
%     % Cheak if x0_i is feasible
%     if isempty(x0_i)
%         continue; 
%     end
% 
%     % Reset configeration
%     config_i = config;
%     config_i.operatingCond.v_avg = v_avg_val;
%     config_i.x0 = x0_i; 
%     config_i.k_h = k_h_val;
%     config_i.k_k = k_k_val;
% 
%     % Generate initial guess trajectory
%     T_i_guess = L_val/v_avg_val;
%     trajInit_i = initial_guess(config, x0_i, T_i_guess);
% 
%     % Run optimization 
%     [Z, P, w, r, constraintNames, h, bounds, boundNames] = getConstraintsAndCosts(config_i, trajInit_i);
%     [solutionIPOPT,  stats] = solveCasadi(config_i, Z, P, r, h, w, bounds, trajInit_i);
% 
%     % Extract results
%     solution = full(solutionIPOPT.x);
%     trajPlot = interpolate_Z2traj(config_i, solution, [], 0, dt);
%     traj_full = trajToFulltraj(config_i, trajPlot);
%     getAnimationRABBIT(traj_full, 0, config, []);
% 
% end



















%% Initialize Simulation Settings
% ode.settings = odeset('RelTol', 1e-7, 'AbsTol', 1e-8);
% ode.solver   = @ode45;
% dt = 3e-2;  % Time step
% 
% % Define parameters to optimize
% optParameterNames = config.optParameterNames;
% 
% %% Test inverse kinematics
% % Bipedal_IK(config, step_length, hip_height, phi, v_target)
% x = Bipedal_IK(config, 0.4, 0.8, 0, 0.4);
% if isempty(x)
%     disp('not feasible')
% end
% 
% %% 1. 定义极其精简的网格
% v_avg_grid = [0.5, 1.0];   % 只测两个速度
% L_grid     = [0.3, 0.4];   % 只测两个步长
% H_grid     = 0.8;          % 固定高度
% phi_grid   = 0;            % 固定躯干
% p_grid     = [50, 100];    % 两个刚度
% 
% % 将网格拍平 (Flattening)
% [V, L, H, P, PHI] = ndgrid(v_avg_grid, L_grid, H_grid, p_grid, phi_grid);
% v_vector = V(:);
% L_vector = L(:);
% H_vector = H(:);
% p_vector = P(:);
% phi_vector = PHI(:);
% 
% num_trials = numel(v_vector); % 总共只有 2*2*1*2*1 = 8 组测试
% 
% %% Grid config
% % Operating Conditions
% v_avg_grid = 0.5:0.1:1.5;
% 
% % Physical Guesses
% L_grid   = 0.1:0.2:0.6;
% H_grid   = 0.75:0.05:0.85;
% phi_grid = -0.05:0.05:0.05;
% 
% % Model parameter grid
% p_grid = logspace(log10(1), log10(150), 5);
% 
% input_grids = {v_avg_grid, L_grid, H_grid, p_grid, phi_grid};
% 
% %% Parallel
% parfor i = 1:numel(v_vector)
%     % Get the parameters of the current grid point
%     v_avg_val = v_vector(i);
%     L_val    = L_vector(i);
%     H_val    = H_vector(i);
%     k_val    = p_vector(i);
%     phi_val  = phi_vector(i);
% 
%     X0_ik = Bipedal_IK(config, L_val, H_val, phi_val, v_avg_val);
% 
%     if isempty(X0_ik)
%         continue; 
%     end
% 
%     config_i = config;
%     config_i.x0 = X0_ik; 
%     config_i.paramValues.k_l = k_val;
% 
%     trajSim = flow(config_i, dt);
% 
% end
% %%
% % Save original configuration and set backup parameters
% config_0 = config;
% backup_batchSize = 1000;
% 
% % Set optimization parameter bounds and simulation values
% config.optParameterUpperBound = [200];
% config.simValue.alpha_k_p = 20;
% config.simValue.alpha_k_d = 5;
% config.simValue.alpha_tau_max = 30;
% 
% % Set hip compliance parameter
% config.paramValues.k_h  = 5;
% 
% %% Initialize Model Parameters and State
% % Get model parameters and initial state
% p_vec = getModelParameters(config, 0, 0);
% X0    = config.x0;
% 
% %% Initialize Simulation Settings
% ode.settings = odeset('RelTol', 1e-7, 'AbsTol', 1e-8);
% ode.solver   = @ode45;
% dt = 3e-2;  % Time step
% 
% % Define parameters to optimize
% optParameterNames = config.optParameterNames;
% 
% %% Define Parameter Grids for Warm Start
% % Velocity grid for operating conditions
% v_avg_grid   = 0.1:0.05:1.5;
% num_v_avg    = numel(v_avg_grid);
% 
% % Define parameter grids for warm start exploration
% p_grid       = logspace(log10(20), log10(150), 6);      % Stiffness
% alpha0_grid  = linspace(0.25, 0.65, 8);                % Initial angle
% dx0_grid     = linspace(0.8, 3, 8);                    % Initial x-velocity
% dy0_grid     = 0;                                      % Initial y-velocity
% u_alpha_grid = -6:3:3;                                 % Angle control input
% u_l_grid     = 3:3:12;                                 % Length control input
% 
% % Create multidimensional grid
% input_grids = {p_grid, u_alpha_grid, u_l_grid, alpha0_grid, dx0_grid, dy0_grid};
% n_g = numel(input_grids);
% grid_outputs = cell(1, n_g);
% [grid_outputs{:}] = ndgrid(input_grids{:});
% 
% % Flatten grids into vectors
% p_vector       = grid_outputs{1,1}(:).';
% u_alpha_vector = grid_outputs{1,2}(:).';
% u_l_vector     = grid_outputs{1,3}(:).';
% alpha0_vector  = grid_outputs{1,4}(:).';
% dx0_vector     = grid_outputs{1,5}(:).';
% dy0_vector     = grid_outputs{1,6}(:).';
% 
% clear grid_outputs;
% 
% %% Initialize Data Storage
% % Test trajectory to determine state dimension
% trajTest = flow(config, dt);
% for j = 1:numel(trajTest)
%     trajTest(j).p = [];
% end
% config.optParameterNames = {};
% Z_test   = interpolate_traj2Z(switch_config(config), trajTest);
% n_z = numel(Z_test);
% 
% % Initialize result arrays
% operatingPoint_array = v_avg_grid';
% cost_array = inf(num_v_avg, 1);
% Z_array = zeros(num_v_avg, n_z);
% warmStart_array = nan(num_v_avg, n_g);
% gridPoint_array = nan(num_v_avg, 2);
% numConverged_array = zeros(num_v_avg, 1);
% largestNorm     = zeros(num_v_avg, 1);
% numberTrials    = zeros(num_v_avg, 1);
% 
% % Initialize tracking arrays
% trackVavg = nan(2, numel(dy0_vector));
% trackCost = nan(2, numel(dy0_vector));
% trackOptParam = nan(2, numel(dy0_vector), size(p_grid, 1));
% trackVavgSim = nan(1, numel(dy0_vector));
% r_f = p_vec(getParameterIndex('r_f', config));
% 
% %% Restart the parallel processing pool. 
% disp('starting warm start')
% delete(gcp('nocreate'));  % close any existing pool
% parpool('local');         % start a fresh pool
% %% Main Optimization Loop
% kk = 0;
% failedIterations = [];
% 
% while kk < numel(dy0_vector)
%     % Preallocate results structure for current batch
%     results = repmat(struct( ...
%         'idx', [], ...
%         'success', false, ...
%         'cost', inf, ...
%         'Z', [], ...
%         'operatingPoint', nan, ...
%         'gridPoint', nan(1,2), ...
%         'warmStart', nan(1,4), ...
%         'trackCost', nan(2,1), ...
%         'trackOptParam', nan(2, size(p_grid, 1)), ...
%         'trackVavg', nan, ...
%         'trackVavgSim', nan, ...
%         'normToOptimal', 0, ...
%         'failReason', "" ...
%     ), 1, min(kk+backup_batchSize, numel(dy0_vector))-kk);
% 
%     % Parallel processing of batch
%     parfor i = 1:min(kk+backup_batchSize, numel(dy0_vector))-kk
%         try 
%             localResult = results(i);  % Local result structure
% 
%             % Configure current iteration
%             config_i = config;
%             config_i.paramValues.k_l = p_vector(i+kk);
%             config_i.optParameterNames = {};
%             p_vec = getModelParameters(config_i, 0, 0);
% 
%             % Set initial conditions and control inputs
%             config_i.x0  = update_initial_cond(config_i, X0, p_vec, alpha0_vector(i+kk), dx0_vector(i+kk), dy0_vector(i+kk));
%             config_i.simValue.constantInput = [u_alpha_vector(i+kk), u_l_vector(i+kk)];
%             config_i_sim = switch_config(config_i);
% 
%             % Simulate trajectory
%             trajSim_FSFS = flow(config_i_sim, dt);
%             trajSim      = switch_traj(trajSim_FSFS);
%             continueFlag = false;
% 
%             % Check simulation validity
%             for phaseIdx = 1:numel(trajSim)
%                 if length(trajSim(phaseIdx).t) < 5 || length(trajSim(phaseIdx).t) > 200 || any(trajSim(phaseIdx).lambda < -1e-2)
%                     continueFlag = true;
%                     localResult.failReason = ['invalid ', config_i.phaseSequence{phaseIdx}, ' in sim'];
%                     break;
%                 end
%             end
% 
%             if continueFlag
%                 results(i) = localResult;
%                 continue;
%             end
% 
%             % Zero out control inputs for interpolation
%             for dd = 1:numel(trajSim)
%                 trajSim(dd).u = 0*trajSim(dd).u;
%             end
%             Z_init = interpolate_traj2Z(config_i, trajSim);
% 
%             % Estimate average velocity from simulation
%             d0    = trajSim(1).x(end,1) + trajSim(1).x(end,3)*r_f  + trajSim(1).x(end,4)* sin(trajSim(1).x(end,3));
%             x_rec = d0 - r_f * (trajSim(2).x(end,1)) - trajSim(2).x(end,2)*sin(trajSim(2).x(end,1));
%             totalTime = sum(arrayfun(@(tr) tr.t(end), trajSim));
%             v_avgSim = (x_rec - trajSim(1).x(1,1))/totalTime;
% 
%             if v_avgSim < 0
%                 localResult.failReason = 'backward motion in sim';
%                 results(i) = localResult;
%                 continue;
%             end
% 
%             % Find closest target velocity
%             localResult.trackVavgSim = v_avgSim;
%             [~, idx] = min(abs(v_avg_grid - v_avgSim));
%             localResult.idx = idx;
%             localResult.trackVavg = v_avg_grid(idx);
%             localResult.operatingPoint = v_avg_grid(idx);
% 
%             % Step 1: Optimize with fixed stiffness
%             config_i.operatingCond.v_avg  = v_avg_grid(idx);
%             [Z_fixed, P_fixed, w_fixed, r_fixed, ~, h_fixed, bounds_fixed, ~] = getConstraintsAndCosts(config_i, trajSim);
% 
%             [solutionIPOPT, stats] = solveCasadi(config_i, Z_fixed, P_fixed, r_fixed, h_fixed, w_fixed, bounds_fixed, trajSim);
% 
%             if ~strcmp(stats.return_status, 'Solve_Succeeded')
%                 failedIterations = [failedIterations; {i, 'opt 1 failed'}];
%                 results(i) = localResult;
%                 continue;
%             end
% 
%             solution = full(solutionIPOPT.x);
%             localResult.trackCost(1) = full(solutionIPOPT.f);
%             localResult.trackOptParam(1,:) = p_vector(i+kk);
% 
%             % Interpolate solution
%             trajSol = interpolate_Z2traj(config_i, solution);
% 
%             % Step 2: Optimize with adaptive stiffness
%             p_matrix = p_vector(i+kk)*ones(size(trajSol(1).x,1), size(p_vector(i+kk),2));
%             [trajSol.p] = deal(p_matrix);
% 
%             config_i.optParameterNames = optParameterNames;
%             config_i.optParameterInit  = p_vector(i+kk);
%             [Z_adaptive, P_adaptive, w_adaptive, r_adaptive, ~, h_adaptive, bounds_adaptive, names] = getConstraintsAndCosts(config_i, trajSol);
% 
%             [solutionIPOPT, stats] = solveCasadi(config_i, Z_adaptive, P_adaptive, r_adaptive, h_adaptive, w_adaptive, bounds_adaptive, trajSol);
% 
%             if ~strcmp(stats.return_status, 'Solve_Succeeded')
%                 localResult.failReason = 'opt 2 failed';
%                 results(i) = localResult;
%                 continue;
%             end
% 
%             solution = full(solutionIPOPT.x);
%             traj     = interpolate_Z2traj(config_i, solution);
%             localResult.trackCost(2) = full(solutionIPOPT.f);
%             localResult.trackOptParam(2,:) = traj(1).p(1,:);
% 
%             % Check if fixed stiffness performed better
%             if localResult.trackCost(1) < localResult.trackCost(2)
%                  localResult.failReason = 'fixed stiffness has better cost';
%                  results(i) = localResult;
%                  continue;
%             end
% 
%             % Store successful result
%             [traj.p] = deal([]);
%             localResult.Z             = interpolate_traj2Z(config_i, traj);
%             localResult.cost          = localResult.trackCost(2);
%             localResult.gridPoint     = [v_avg_grid(idx), localResult.trackOptParam(2, :)];
%             localResult.warmStart     = [alpha0_vector(i+kk), dx0_vector(i+kk), dy0_vector(i+kk), u_alpha_vector(i+kk), u_l_vector(i+kk), p_vector(i+kk)];
%             localResult.normToOptimal = norm(Z_init-localResult.Z)^2;
%             localResult.success       = true;
% 
%             results(i) = localResult;
% 
%         catch ME
%             localResult.failReason = ['unhandled error: ', ME.identifier];
%             results(i) = localResult;
%             continue;
%         end
%     end
% 
%  % Process results from current batch
%     for i = 1:numel(results)
%         if ~results(i).success || numel(results(i).failReason)>1
%             failedIterations = [failedIterations; {i+kk, results(i).failReason}];
%             continue;
%         elseif isempty(results(i).idx)
%             failedIterations = [failedIterations; {i+kk, 'no idx found'}];
%             continue;
%         end
% 
%         idx = results(i).idx;
%         numberTrials(idx) = numberTrials(idx) + 1;
% 
%         % Update best solution if current is better
%         if results(i).success && results(i).cost < cost_array(idx)
%             Z_array(idx,:)         = results(i).Z;
%             cost_array(idx)        = results(i).cost;
%             gridPoint_array(idx,:) = results(i).gridPoint;
%             warmStart_array(idx,:) = results(i).warmStart;
%             operatingPoint_array(idx,:) = results(i).operatingPoint;
%             numConverged_array(idx) = 1;
%             largestNorm(idx) = results(i).normToOptimal;
%         % Handle similar solutions
%         elseif norm(results(i).cost - cost_array(idx)) < 1e-2 && norm(Z_array(idx,:) - results(i).Z') < 1e-1 && norm(gridPoint_array(idx,2)-results(i).gridPoint(2)) < 1   
%             numConverged_array(idx) = numConverged_array(idx)+ 1;
%             if results(i).normToOptimal > largestNorm(idx)
%                 largestNorm(idx) = results(i).normToOptimal;
%                 warmStart_array(idx,:) = results(i).warmStart;
%             end
%         end
% 
%         % Update tracking arrays
%         trackCost(:, i+kk)        = results(i).trackCost;
%         trackOptParam(:, i+kk, :) = results(i).trackOptParam;
%         trackVavg(i+kk)           = results(i).trackVavg;
%         trackVavgSim(i+kk)        = results(i).trackVavgSim;
%     end 
% 
%     % Save intermediate results
%     save_data([dated_filename, '_UNFINISHED'], 'o', 'Z_array', Z_array, 'gridPoint_array', gridPoint_array, 'cost_array', cost_array, ...
%         'warmStart_array', warmStart_array, 'numConverged_array', numConverged_array ,'largestNorm', largestNorm , 'numberTrials', numberTrials, 'config', config);
% 
%     % Move to next batch
%     kk = kk+backup_batchSize;
% end
% %% Save Final Results and Analysis
% diagnosis = {};
% diagnosis.alpha0_grid      = alpha0_grid;
% diagnosis.p_grid           = p_grid;
% diagnosis.dx0_grid         = dx0_grid;
% diagnosis.dy0_grid         = dy0_grid;
% diagnosis.failedIterations = failedIterations;
% diagnosis.trackCost        = trackCost;
% diagnosis.trackOptParam    = trackOptParam;
% diagnosis.trackVavg        = trackVavg;
% 
% % Save final data
% filename = save_data(dated_filename, 'v', 'Z_array', Z_array, 'gridPoint_array', gridPoint_array, 'cost_array', cost_array, ...
%     'warmStart_array', warmStart_array, 'numConverged_array', numConverged_array ,'largestNorm', largestNorm , 'diagnosis', diagnosis, 'numberTrials', numberTrials, 'config', config_0);
% 
% % Analyze and document results
% analyzeFailuresWarmStart(failedIterations, numel(dy0_vector), [filename,'.txt']);
% documentInequalityConstraints(config, [filename,'_ineq.txt'])
% 
% disp('warm start finished')
% 
% %% Visualization and Analysis Section
% % Example: Plot solution for specific velocity index
% idx_toPlot = 31;
% config.optParameterNames = {};
% config.paramValues.k_l = gridPoint_array(idx_toPlot,2);
% trajSol   = interpolate_Z2traj(config, Z_array(idx_toPlot,:)', [], 1, 3e-2);
% 
% % Create animation and plots
% traj_full = map2ClassTrajectory(config, trajSol);
% frames    = getAnimationHopper5DOF(traj_full,1,false);
% 
% % Plot trajectory with constraints
% [activeBounds, trajLB, trajUB] = checkActiveIneqConstraints(trajSol, config);
% plot_traj(config, trajSol, [],  trajLB, trajUB);
% 
% % Verify average velocity
% v_avg_test = traj_full(end).x(end,1)/(traj_full(1).t(end)+traj_full(2).t(end));
% disp(v_avg_test);
% 
% % Analyze optimization error
% error = analyze_error(config, Z_array(idx_toPlot,:)', 1e-4, 1);
% 
% %% Helper Functions
% function X0_updated = update_initial_cond(config, X0, p_vec, alpha0, dx0, dy0)
%     % Update initial conditions for simulation
%     % Inputs:
%     %   config - system configuration
%     %   X0 - base initial state vector
%     %   p_vec - parameter vector
%     %   alpha0 - initial angle
%     %   dx0 - initial x-velocity
%     %   dy0 - initial y-velocity
% 
%     l_0 = p_vec(getParameterIndex('l_0', config));
%     X0_updated = X0;  % has 8 elements
% 
%     X0_updated(1) = 0;        % x-position
%     X0_updated(2) = 1.5;      % y-position (just enough clearance)
%     X0_updated(3) = alpha0;   % angle
%     X0_updated(4) = l_0;      % leg length
% 
%     X0_updated(5) = dx0;      % x-velocity
%     X0_updated(6) = dy0;      % y-velocity
%     X0_updated(7) = 0;        % angular velocity
%     X0_updated(8) = 0;        % leg velocity
% end
% 
% function traj_min = switch_traj(traj)
%     % Switch trajectory phases and adjust coordinates
%     % Assumes trajectory is generated using HF_noPitch with sequence Stance then Flight
%     % Output: modified trajectory with switched phases and reduced states
% 
%     traj_min = traj(end-1:end);
%     traj_min(1).x(:,1) = traj_min(1).x(:,1) - traj_min(1).x(1,1);
% end
% 
% function config_HF = switch_config(config)
%     % Switch configuration for different phase sequence
%     config_HF = config;
%     config_HF.phaseSequence = {'Flight', 'Stance', 'Flight', 'Stance'};
%     config_HF.N             = [config.N, config.N];
%     config_HF.simControl    = {'PD', 'constant', 'PD', 'zero'};
%     config_HF.simValue.alpha_target = config_HF.x0(3);
% end
% 
% 
% 
% 
% 






