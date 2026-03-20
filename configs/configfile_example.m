% Example Configuration for hopping optimization with adaptive hip actuation
% Sets up model parameters, optimization constraints, and continuation settings

config = struct();
%% Model Configuration
config.model_name = 'HF_compHip_min';          % Composite hip actuation model
config.phaseSequence = {'Flight','Stance'};    % Phase sequence for hopping gait
config.N = [50, 50];                           % Collocation segments per phase

%% Operating Conditions
config.operatingCond.v_avg = 0.6;              % Average forward velocity [m/s]
config.operatingCond.x0    = 0;                % Initial x-position [m]

%% Numerical Discretization
config.collocationScheme = 'hermiteSimpson';   % Collocation method for dynamics
config.inputInterpolation = 'piecewiseLinear'; % Control input interpolation

%% Periodicity Constraints
config.periodicityFlag = 1;
% Matrix defining which states are periodic (excludes x-position for forward hopping)
config.periodicityMatrix =   [0 1 0 0 0 0 0 0;...     % y-position
                              0 0 1 0 0 0 0 0;...     % alpha (torso angle)  
                              0 0 0 1 0 0 0 0;...     % l (leg length)
                              0 0 0 0 1 0 0 0;...     % dx
                              0 0 0 0 0 1 0 0;...     % dy
                              0 0 0 0 0 0 1 0;...     % dalpha
                              0 0 0 0 0 0 0 1];       % dl

%% Optimization Parameters
config.optParameterNames = {'k_l'};                    % Leg stiffness to optimize
config.optParameterInit  = [20];                       % Initial guess
config.optParameterLowerBound = [1];                   % Minimum leg stiffness
config.optParameterUpperBound = [200];                 % Maximum leg stiffness

% Alternative: No parameter optimization
% config.optParameterNames = {};
% config.optParameterInit  = [];
% config.optParameterLowerBound = [];
% config.optParameterUpperBound = [];

%% Optimize time
config.optimizeTimeFlag = 1;    % for now always 1

%% Dynamics and Constraints Configuration
config.GRFfunc = '.Stance.S2F.event';              % Define function in the model foe Ground reaction force calculation 
config.useInequalityConstraints = true;            % Enable inequality constraints

%% Cost Function Configuration
config.costName = 'weighted_u_squared';            % Type of cost function to minimize
config.costNormalization = 'CoT_v_avg';            % Cost normalization method

%% Initial Guess Generation
config.loadInitialGuessFlag = 0;                   % Flag to load initial guess from file
config.loadInitialGuessName = '';                  % Filename for initial guess data

%% Forward Simulation Settings for Initial Guess 
% these settings are used in the flow map and are used specififc
config.simInitialGuessFlag = 1;                    % Generate initial guess via forward simulation
config.simControl = {'PD', 'constant'};            % Control types for each phase
config.simValue.constantInput = [1, -1];           % Constant control input values
config.simValue.alpha_target = 0.3500;             % Target hip angle for PD control

% Initial state vector [x, y, alpha, l, dx, dy,  dalpha, dl]
config.x0 = [0; 0.855723557545085; -0.276089751432568; 0.837438435350479; ...
             1.449634584498169; 0.851354792204313; -1.312931214136162; 1.196381903510173];    

%% Variable Mapping and Default Values to class variables
config.defaultStateValues = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];  % Default state values
config.defaultInputValues = [0; 0];                           % Default input values
config.defaultStateNames  = {'x', 'y', 'phi', 'alpha', 'l', 'dx', 'dy','dphi', 'dalpha', 'dl'};
config.defaultInputNames  = {'u_alpha', 'u_l'};

%% Post-processing and Output Options
% for now unused
config.saveResultsFlag = 0;                        % Save optimization results flag
config.plotResultsFlag = 0;                        % Plot results flag
config.startContinuationFlag = 0;                  % Start parameter continuation flag

%% Continuation Study Configuration
config.cont.type = 'Poorman';                      % Continuation method type

% Operating condition grid settings
config.cont.gridOpCondNames = {'v_avg'};           % Operating condition variables
config.cont.gridOpCondMin   = [0.1];               % Minimum operating condition values
config.cont.gridOpCondMax   = [3];                 % Maximum operating condition values
config.cont.gridOpCondStep  = [0.025];             % Operating condition step sizes

% Parameter grid settings
config.cont.gridParamNames = {'k_l'};              % Parameter variables for continuation
config.cont.gridParamMin   = [1];                  % Minimum parameter values
config.cont.gridParamMax   = [200];                % Maximum parameter values
config.cont.gridParamStep  = [1];                  % Parameter step sizes

% Cost acceptance thresholds
config.cont.costSlack = 1e-5;                      % Cost improvement threshold for replacement
config.cont.costNew = 1e-3;                        % Cost threshold for new solution acceptance

% Batch processing settings
config.cont.n_batch = 500;                         % Number of points per batch
config.cont.n_out   = 5*config.cont.n_batch;       % Maximum output solutions

% File management for continuation
config.cont.filenameMAT = 'PC';                    % Base filename for continuation data
config.cont.warmStartMAT = '21_08_WS_weighted_u_squared_CoT_v_avg_HF_compHip_min_kh10';

% Visualization settings for continuation results
config.cont.plotViews = {[1,2]};                   % Parameter pairs to visualize
config.cont.sliceValues = [];                      % Fixed values for slicing high-dimensional data
config.cont.sliceTolerance = 1e-5;                 % Tolerance for slice matching
