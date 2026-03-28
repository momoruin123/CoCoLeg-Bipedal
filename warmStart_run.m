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
% configfile_BF_min;
configfile_BF_run;

config.optConfig.print_level = 0;

% Save original configuration and set backup parameters
config_0 = config;
backup_batchSize = 50;

%% Initialize Model Parameters and State
% Get model parameters and initial state
p_vec = getModelParameters(config, 0, 0);
x0    = config.x0;

% Set hip compliance parameter
config.paramValues.k_h = 20;
config.paramValues.k_k = 20;
%% Initialize Simulation Settings
ode.settings = odeset('RelTol', 1e-7, 'AbsTol', 1e-8);
ode.solver   = @ode45;
dt = 3e-2;  % Time step

% Define parameters to optimize
optParameterNames = config.optParameterNames;

%% Define Parameter Grids for Warm Start
% grid sample of velocity
% v_avg_grid   = [1.4, 1.7];
% v_avg_grid   = [0.01, 0.1:0.1:1.5];
v_avg_grid   =0.5:0.5:5;

% v_avg_grid   = 0.001:0.1:1.501;
num_v_avg    = numel(v_avg_grid);

% Latin hypercube sample of params
% k_h_grid, k_k_grid, L_grid, H_grid, phi_grid

% Config. of sampling
n_samples = 20;
% lb = [1e-6, 1e-6, 0, 0, -pi/4]; % lower bounds
% ub = [1000, 1000, 1.2, 0.7, pi/4]; % upper bounds

% k_h_grid, k_k_grid, T, percentage
lb = [1e-6, 1e-6, 0.1, 0.4]; % lower bounds
ub = [300, 300, 0.6, 0.8]; % upper bounds

n_g = numel(lb);
n_p = numel(optParameterNames);    % number of optimized params

% Latin_hypercube_sampling
input_matrix = latin_hypercube_sampling(lb, ub, n_samples);
k_h_grid   = input_matrix(:,1);     % hip stiffness
k_k_grid   = input_matrix(:,2);     % knee stiffness
T_grid     = input_matrix(:,3);
percentage_grid = input_matrix(:,4);

% L_grid     = input_matrix(:,3);     % step length
% H_grid     = input_matrix(:,4);     % hip height
% phi_grid   = input_matrix(:,5);     % pitch angle

% Flatten grids into vectors
input_matrix = repmat(input_matrix, num_v_avg, 1);

% Define parameter grids for warm start exploration
k_h_vector   = input_matrix(:,1);     % hip stiffness
k_k_vector   = input_matrix(:,2);     % knee stiffness
T_vector     = input_matrix(:,3);
percentage_vector   = input_matrix(:,4);

% L_vector     = input_matrix(:,3);     % step length
% H_vector     = input_matrix(:,4);     % hip height
% phi_vector   = input_matrix(:,5);     % pitch angle
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
trackOptParam = nan(2, n_p, numel(v_avg_vector)); 

%% Restart the parallel processing pool. 
disp('starting warm start')
delete(gcp('nocreate'));  % close any existing pool
setenv('OMP_NUM_THREADS', '1');
setenv('MKL_NUM_THREADS', '1');
parpool('local');  % Number of local core

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
    parfor i = 1:min(kk+backup_batchSize, numel(v_avg_vector))-kk
        localResult = results(i);

        try
            % Configure current iteration
            config_i = config;
            config_i.paramValues.k_h = k_h_vector(i+kk);
            config_i.paramValues.k_k = k_k_vector(i+kk);
            config_i.optParameterNames = {};
            p_vec = getModelParameters(config_i, 0, 0);

            % Set initial conditions and control inputs
            % config_i.x0 = Bipedal_IK(config, L_vector(i+kk), H_vector(i+kk), phi_vector(i+kk), v_avg_vector(i+kk));
            config_i.x0 = x0;

            % Cheak physical feasibility
            if isempty(config_i.x0)
                localResult.failReason = 'IK infeasible';
                results(i) = localResult;
                continue;
            end

            % Generate initial guess trajectory
            % 1. Interpolation
            T_i_guess = T_vector(i+kk);
            percentage = percentage_vector(i+kk);
            trajInit = initial_guess(config_i, T_i_guess, percentage, 0);
            Z_init = interpolate_traj2Z(config_i, trajInit);

            % 2. Simulation
            % trajInit = flow(config, dt);

            % =========================NEED ADJUST========================
            % Find closest target velocity
            [~, idx] = min(abs(v_avg_grid - v_avg_vector(i+kk)));
            localResult.idx = idx;
            localResult.trackParams = [v_avg_grid(idx); k_h_vector(i+kk); k_k_vector(i+kk); percentage_vector(i+kk); T_vector(i+kk)];
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
            localResult.warmStart     = [k_h_vector(i+kk), k_k_vector(i+kk), percentage_vector(i+kk), T_vector(i+kk)];
            localResult.normToOptimal = norm(Z_init-localResult.Z)^2;
            localResult.success       = true;

        catch ME
            localResult.failReason = ['unhandled error: ', ME.identifier];
            results(i) = localResult;
            continue;
        end
                    
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

%% Save Final Results and Analysis
diagnosis = {};
diagnosis.k_h_grid         = k_h_grid;
diagnosis.k_k_grid         = k_k_grid;
diagnosis.percentage_grid  = percentage_grid;
diagnosis.T_grid           = T_grid;
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
