% config file. Choose all settings for the optimization problem.
config = struct();
%% Choose your model
config.model_name = 'HP_min';

%% Choose you sequence
config.phaseSequence = {'Flight', 'Stance', 'Flight'};
config.N = [30, 30, 30];                    % define the number of segments in each phase for the collocation.


%% Define operating conditions ()
config.operatingCond.HP = 1.5;
config.operatingCond.idxHP = 1;

%% Exclude elements from periodicity
config.periodicityFlag = 1;
% Matrix defining which states are periodic (excludes x-position for forward hopping)
config.periodicityMatrix = eye(4);

%% Define parameters to be optimized on
config.optParameterNames = {'k_l'};
config.optParameterInit  = [20];
config.optParameterLowerBound = [1];
config.optParameterUpperBound = [200];

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
% Define bounds to be used also when the parameter is not optimized upon.
% this is needed for consistent inequualities in a continuation.
config.inequalityParamBoundName = {'k_l'};
config.inequalityParamLowerBound = [1];
config.inequalityParamUpperBound = [200];

%% Define cost and normalization factor:
config.costName = 'u_squared';              % 'u_squared', 'none', 'positive_mechanical_work'
config.costNormalization = 'totalTime';    %  'CoT_v_avg', 'totalTime', 'None';

%% Define initial guess
config.loadInitialGuessFlag = 0;                % if true, load inital guess fro file

%% Forward Simulation Settings for Initial Guess 
% these settings are used in the flow map and are used specififc
config.loadInitialGuessName = '';
% if it is a forward simulation
config.simInitialGuessFlag = 1;                 % if true, run forward simulation to get initial guess
config.simControl = {'zero', 'zero', 'zero'};   % should follow the phase sequence pattern.
config.x0         = [1.5; 1; -1; 0];            % should match the starting phase


%% Variable Mapping and Default Values to class variables
config.defaultStateValues = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
config.defaultInputValues = [0; 0];

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


