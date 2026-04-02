% Example Configuration for bipedal robot optimization with adaptive hip actuation
% Sets up model parameters, optimization constraints, and continuation settings
config = struct();
config.dt = 3e-2;

%% Optimization Configuration
config.optConfig.tolerance = 1;
config.optConfig.print_level = 5;

%% Model Configuration
config.model_name = 'BF_min';          % Composite hip actuation model
config.phaseSequence = {'SingleStance'};    % Phase sequence for hopping gait
% config.N = [50, 50, 50];                           % Collocation segments per phase
config.N = 50;

%% Operating Conditions
config.operatingCond.x0    = 0;                % Initial x-position [m]
% Geit optimizing methods
% config.optGaitMethods = 'Fix_step_length'; 
% config.stepLength = 0.5;
config.optGaitMethods = 'Fix_average_velocity';
config.operatingCond.v_avg = 0.6;             % Average forward velocity [m/s]

%% Numerical Discretization
config.collocationScheme  = 'hermiteSimpson';   % Collocation method for dynamics
config.inputInterpolation = 'piecewiseLinear'; % Control input interpolation

%% Periodicity Constraints
config.periodicityFlag = 1;

% Matrix defining which states are periodic (excludes x-position for forward hopping)
config.periodicityMatrix = eye(10);
% to be added config.periodicityFlag  to signiize if we want to add periodicity constraInt.

%% Optimization Parameters
% config.optParameterNames = {'alpha0', 'beta0', 'k_h', 'k_k'};               % Leg stiffness to optimize
% config.optParameterInit  = [0, 0, 20, 20];                          % Initial guess
% config.optParameterLowerBound = [-pi/2, -pi/2, 1e-6, 1e-6];               % Minimum leg stiffness
% config.optParameterUpperBound = [pi/2, pi/2, 1000, 1000];                 % Maximum leg stiffness

config.optParameterNames = {'k_h', 'k_k'};               % Leg stiffness to optimize
config.optParameterInit  = [20, 20];                          % Initial guess
config.optParameterLowerBound = [1e-6, 1e-6];               % Minimum leg stiffness
config.optParameterUpperBound = [1000, 1000];   

% config.optParameterNames = {'k_h'};               % Leg stiffness to optimize
% config.optParameterInit  = [20];                          % Initial guess
% config.optParameterLowerBound = [1e-6];               % Minimum leg stiffness
% config.optParameterUpperBound = [1000];                 % Maximum leg stiffness

% Alternative: No parameter optimization
% config.optParameterNames = {};
% config.optParameterInit  = [];
% config.optParameterLowerBound = [];
% config.optParameterUpperBound = [];

%% Optimize time
config.optimizeTimeFlag = 1;    % for now always 1
% config.Time = 0.7;              % if don't optimize time, set time length

%% Dynamics and Constraints Configuration
% config.GRFfunc = '.SingleStance.S2S.jump';          % Define function in the model foe Ground reaction force calculation 
config.GRFfunc = '.SingleStance.compute_lambda';          % Define function in the model foe Ground reaction force calculation 
config.useInequalityConstraints = true;            % Enable inequality constraints

% Define bounds to be used also when the parameter is not optimized upon.
% this is needed for consistent inequualities in a continuation.
% config.inequalityParamBoundName = {'k_h'};
% config.inequalityParamLowerBound = [1];
% config.inequalityParamUpperBound = [200]; 

q_lb = [-0.2*pi; -0.5*pi; -0.5*pi; -pi/2; -pi/2];
q_ub = [0.2*pi;  0.5*pi;  0.5*pi; -0.001; -0.001];
dq_lb = [-12; -12; -12; -15; -15];
dq_ub = [12; 12; 12; 15; 15];

% ub_state = [0.5*pi;  inf;  inf; inf; inf];
% lb_state = [-0.5*pi; -inf; -inf; -inf; -inf];
config.q_lb  = q_lb;
config.q_ub  = q_ub;
config.dq_lb = dq_lb;
config.dq_ub = dq_ub;

config.input_lb = [-300; -300; -220; -220];
config.input_ub = [300; 300; 220; 220];

config.dt_bounds = [0.002, 0.05];

%% Cost Function Configuration
% none; u_squared; weighted_u_squared; positive_mechanical_work
config.costName = 'u_squared';            % Type of cost function to minimize
config.costNormalization = 'CoT_v_avg';            % Cost normalization method

%% Initial Guess Generation
config.loadInitialGuessFlag = 0;                   % Flag to load initial guess from file
config.loadInitialGuessName = '';                  % Filename for initial guess data

%% Forward Simulation Settings for Initial Guess 
% these settings are used in the flow map and are used specififc
config.simInitialGuessFlag = 1;                    % Generate initial guess via forward simulation
config.simControl = {'zero'};                      % Control types for each phase
config.simValue.constantInput = [2, -1, 1, -0.5];  % Constant control input values
config.simValue.alpha_target = 0.3500;             % Target hip angle for PD control
config.simValue.alpha_k_p = 2.5;
config.simValue.alpha_k_d = 1;
config.simValue.alpha_tau_max = 30;

% Initial state vector.
% Swing leg landing initial codiation
% config.x0 = [-0.1; 0.29; 0.21; -0.1; -0.54; ...
%              0;    -3; 3;    0;    -1.5];
% Swing leg landing initial codiation
% config.x0 = [-0.1; 0.29; -0.21; -0.1; -0.34; ...
%              0;    0; 0;    0;    0]; 
% config.x0 = [-0.1; 0.29; 0.21; -0.1; -0.54; ...
%              0;    -3; 3;    0;    -1.5]; 
% config.x0 = [-0.1; 0.29; 0.21; -0.1; -0.54; ...
%              0;    -0; 0;    0;    0];  
% config.x0 = [-pi*0.5; 0.29; 0.21; -0.1; -0.54; ...
%              0;    -3; 3;    0;    -1.5];  
% config.x0 = [0; -0.06; 0.4; -0.31; -0.4; ...
%              0;  -1; 0.5;    0;    0]; 
% config.x0 = [0;0;0; -0.06; 0.4; -0.31; -0.4; ...
%              0;0;0;  -3; 0.1;    0;    0];     
% config.x0 = [0;0;0; -0.06; 0.4; -0.31; -0.5; ...
%              0;0;0;  -3; -1;    0;    0];     

% Walking forward initial codiation (from DCLL)
% config.x0 = [-0.1; 0.29; 0.21; -0.1; -0.54; ...
%              0;    -0.2; 0.65;    0;    -0.9];   

% For optimization
config.x0 = [0; 0.7;-0.5;-0.4;-0.1;
             0;  0;   0; 0;  0];
%% Variable Mapping and Default Values to class variables
config.defaultStateValues = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];  % Default state values
config.defaultInputValues = [0; 0; 0; 0];                    % Default input values
config.defaultStateNames  = {'phi', 'alpha_l', 'alpha_r', 'beta_l', 'beta_r', ...
                             'dphi', 'dalpha_l', 'dalpha_r', 'dbeta_l', 'dbeta_r'};
config.defaultInputNames  = {'u_hl', 'u_hr', 'u_kl', 'u_kr'};

%% Post-processing and Output Options
% for now unused
config.saveResultsFlag = 0;                        % Save optimization results flag
config.plotResultsFlag = 0;                        % Plot results flag
config.startContinuationFlag = 0;                  % Start parameter continuation flag

%% Continuation Study Configuration
config.cont.type = 'Poorman';                      % Continuation method type

% Operating condition grid settings
config.cont.gridOpCondNames = {};           % Operating condition variables
config.cont.gridOpCondMin   = [0.1];               % Minimum operating condition values
config.cont.gridOpCondMax   = [3];                 % Maximum operating condition values
config.cont.gridOpCondStep  = [0.025];             % Operating condition step sizes

% Parameter grid settings
config.cont.gridParamNames = {'k_h', 'k_k'};              % Parameter variables for continuation
config.cont.gridParamMin   = [0, 0];                  % Minimum parameter values
config.cont.gridParamMax   = [20, 20];                % Maximum parameter values
config.cont.gridParamStep  = [1, 1];                  % Parameter step sizes

% Cost acceptance thresholds
config.cont.costSlack = 1e-5;                      % Cost improvement threshold for replacement
config.cont.costNew = 1e-3;                        % Cost threshold for new solution acceptance

% Batch processing settings
config.cont.n_batch = 500;                         % Number of points per batch
config.cont.n_out   = 5*config.cont.n_batch;       % Maximum output solutions

% File management for continuation
config.cont.filenameMAT = 'PC';                    % Base filename for continuation data
config.cont.warmStartMAT = '13_02_WS_u_squared_CoT_v_avg_BF_min_11';

% Visualization settings for continuation results
config.cont.plotViews = {[1,2]};                   % Parameter pairs to visualize
config.cont.sliceValues = [];                      % Fixed values for slicing high-dimensional data
config.cont.sliceTolerance = 1e-5;                 % Tolerance for slice matching