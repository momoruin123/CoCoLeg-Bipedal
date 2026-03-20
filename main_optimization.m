% main file for running an optimization.
% this file will first run a forward simulation of the hopper and then use it to initilize the optimization.
% The file examplifies this in three exmaples: HP_min, HF_full, HF_noPitch_min 
%
%   author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025
%% HP_min: Hopping in place with minimum coordiantes
configfile_HP;
dt = 3e-2;
traj = flow(config, dt);
% Animate forward simulation
trajAnim = map2ClassTrajectory(config, traj);
frames = getAnimationHopper5DOF(trajAnim,1,false);
plot_traj(config, traj);

% Assign the operating cond to the same mean hopping height
config.operatingCond.HP = (mean(traj(1).x(:,1)) + mean(traj(3).x(:,1)))/2;
% Set up and solve optimization
[Z, P, w, r, constraintNames, h, bounds, boundNames] = getConstraintsAndCosts(config, traj);
[solutionIPOPT,  stats] = solveCasadi(config, Z, P, r, h, w, bounds, traj);

% Process solution
solution  = full(solutionIPOPT.x);
trajSol   = interpolate_Z2traj(config, solution, [], 1, dt);
traj_full = map2ClassTrajectory(config, trajSol);
frames    = getAnimationHopper5DOF(traj_full,1,false);

% Analyze results
trajPlot = interpolate_Z2traj(config, solution, [], 0, dt);
[activeBounds, trajLB, trajUB] = checkActiveIneqConstraints(trajPlot, config);
plot_traj(config, trajPlot, [], trajLB, trajUB);
error = analyze_error(config, solution, 1e-4, 1);

% Check optimality conditions
[jacobians_Sym, optimalityConditions_Num] = checkOptimalityConditions(config, solutionIPOPT, Z, P, w, r, h);
const_jacobian          = [optimalityConditions_Num.dr_dZ; optimalityConditions_Num.dh_dZ];
violdatedConstraintName = constraintNames(optimalityConditions_Num.violatedEqualityConstIdx);
activeContraintName     = boundNames(optimalityConditions_Num.activeInequalityConstIdx);

% Visualize sparsity patterns
figure; spy(const_jacobian);
figure; spy(optimalityConditions_Num.ddL_dZ);
disp(['second order optimality (smallest eigenvalue) ', num2str(optimalityConditions_Num.minEigSecondOrder)]);


%% HF_full
% Load config file
configfile_HF_full;
% turn off optimizing the parameter
config.optParameterNames = {};
% Modify parameters
config.paramValues.k_l  = 50;  % change the initial value for the optimized variable
config.paramValues.k_h  = 5;
% Modify the initial conditions to a post impact state
impactFcn = str2func([config.model_name, '.Flight', '.F2S', '.jump']);
P = getModelParameters(config, 0, 0);
config.x0 = impactFcn(config.x0, [], P);
config.x0(2) = config.x0(5)*cos(config.x0(4)) + P(getParameterIndex('r_f', config));
% Simulate system
dt = 3e-2;
traj = flow(config, dt);
frames   = getAnimationHopper5DOF(traj,1,false);

% Calculate average velocity
v_avg_test = traj(end).x(end,1)/(traj(1).t(end)+traj(2).t(end));
disp(['v_avg in simulation : ', num2str(v_avg_test)]);
config.operatingCond.v_avg = round(v_avg_test,1);

% Set up and solve optimization
[Z, P, w, r, constraintNames, h, bounds, boundNames] = getConstraintsAndCosts(config, traj);
[solutionIPOPT,  stats] = solveCasadi(config, Z, P, r, h, w, bounds, traj);

% Process solution
solution  = full(solutionIPOPT.x);
trajSol   = interpolate_Z2traj(config, solution, [], 1, dt);
frames    = getAnimationHopper5DOF(trajSol,1,false);

% Analyze results
trajPlot = interpolate_Z2traj(config, solution, [], 0, dt);
[activeBounds, trajLB, trajUB] = checkActiveIneqConstraints(trajPlot, config);
plot_traj(config, trajPlot, [], trajLB, trajUB);

error = analyze_error(config, solution, 1e-4, 1);

% Display final metrics
v_avg_test = trajSol(end).x(end,1)/(trajSol(1).t(end)+trajSol(2).t(end));
disp(['v_avg of the optimal cond: ', num2str(v_avg_test)]);
disp(['cost: ', num2str(full(solutionIPOPT.f))]);

% Check optimality conditions
[jacobians_Sym, optimalityConditions_Num] = checkOptimalityConditions(config, solutionIPOPT, Z, P, w, r, h);
const_jacobian          = [optimalityConditions_Num.dr_dZ; optimalityConditions_Num.dh_dZ];
violdatedConstraintName = constraintNames(optimalityConditions_Num.violatedEqualityConstIdx);
activeContraintName     = boundNames(optimalityConditions_Num.activeInequalityConstIdx);

% Visualize sparsity patterns
figure; spy(const_jacobian);
figure; spy(optimalityConditions_Num.ddL_dZ);
disp(['second order optimality (smallest eigenvalue) ', num2str(optimalityConditions_Num.minEigSecondOrder)]);


%% HF_noPitch_min
% Load config file
configfile_HF_noPitch_min;
% Modify parameters
config.paramValues.k_l  = 50;
config.paramValues.k_h  = 5;
config.simValue.constantInput = [-1,5];

% Simulate system
config_i_sim = switch_config(config);
dt = 3e-2;
trajSim_FSFS = flow(config_i_sim, dt);
trajAnim = map2ClassTrajectory(config_i_sim, trajSim_FSFS);
frames   = getAnimationHopper5DOF(trajAnim,1,false);

% Calculate average velocity
traj  = switch_traj(trajSim_FSFS);
v_avg_test = trajAnim(end).x(end,1)/(traj(1).t(end)+traj(2).t(end));
disp(['v_avg in simulation : ', num2str(v_avg_test)]);
config.operatingCond.v_avg = round(v_avg_test,1);

% Set up and solve optimization
[Z, P, w, r, constraintNames, h, bounds, boundNames] = getConstraintsAndCosts(config, traj);
[solutionIPOPT,  stats] = solveCasadi(config, Z, P, r, h, w, bounds, traj);

% Process solution
solution  = full(solutionIPOPT.x);
trajSol   = interpolate_Z2traj(config, solution, [], 1, 1e-2);
traj_full = map2ClassTrajectory(config, trajSol);
frames    = getAnimationHopper5DOF(traj_full,1,false);

% Analyze results
trajPlot = interpolate_Z2traj(config, solution, [], 0, 1e-3);
[activeBounds, trajLB, trajUB] = checkActiveIneqConstraints(trajPlot, config);
plot_traj(config, trajPlot, [], trajLB, trajUB);

error = analyze_error(config, solution, 1e-4, 1);

% Display final metrics
v_avg_test = traj_full(end).x(end,1)/(traj_full(1).t(end)+traj_full(2).t(end));
disp(['v_avg of the optimal cond: ', num2str(v_avg_test)]);
disp(['cost: ', num2str(full(solutionIPOPT.f))]);

% Check optimality conditions
[jacobians_Sym, optimalityConditions_Num] = checkOptimalityConditions(config, solutionIPOPT, Z, P, w, r, h);
const_jacobian          = [optimalityConditions_Num.dr_dZ; optimalityConditions_Num.dh_dZ];
violdatedConstraintName = constraintNames(optimalityConditions_Num.violatedEqualityConstIdx);
activeContraintName     = boundNames(optimalityConditions_Num.activeInequalityConstIdx);

% Visualize sparsity patterns
figure; spy(const_jacobian);
figure; spy(optimalityConditions_Num.ddL_dZ);
disp(['second order optimality (smallest eigenvalue) ', num2str(optimalityConditions_Num.minEigSecondOrder)]);

function config_HF = switch_config(config)
    % Switch to flight-stance-flight-stance phase sequence
    config_HF = config;
    config_HF.phaseSequence         = {'Flight', 'Stance', 'Flight', 'Stance'};
    config_HF.N                     = [config.N, config.N];
    config_HF.simControl            = {'PD', 'constantAll', 'PD', 'zero'};
    config_HF.simValue.alpha_target = config_HF.x0(3);
end

function traj_min = switch_traj(traj)
    % Extract last two phases and adjust x-position
    traj_min = traj(end-1:end);
    traj_min(1).x(:,1) = traj_min(1).x(:,1) - traj_min(1).x(1,1);
end
