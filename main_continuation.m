%% Parameter Continuation for Hopping Robot Optimization
%
% This script performs parameter continuation to map optimal solutions across a range of operating conditions and system parameters.
%
% METHODOLOGY:
% The continuation method systematically explores the parameter space by:
% 1. Starting from known good solutions seed points, warm starts)
% 2. Exploring the map in a user defined fashion (poorman, ....)
% 3. Building a complete map of optimal trajectories across the grid
%
% PROCESS FLOW:
% - Initialize from warm start data or existing continuation results
% - Set up parameter grid with specified ranges and step sizes
% - Process points in batches using parallel computing
% - Dynamically update progress with real-time visualization
% - Handle failures and maintain solution quality
% - Save results at regular intervals for crash recovery
%
% KEY FEATURES:
% - Robust handling of optimization failures
% - Solution propagation through parameter space
% - Real-time progress monitoring with multiple visualization views
% - Automatic crash recovery with emergency saving
% - Parallel processing for computational efficiency
%
% OUTPUT: Saved data of complete mapping of optimal trajectories across the specified
% parameter ranges, enabling analysis of solution sensitivity and
% performance trade-offs.
%
% author: C. David Remy, Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

%%
clear; clc;
dated_filename = append_date('PC');
% main continutation file
configfile_BF_min;
% save initial config
config_0 = config;
% change some settings in the config
config.optParameterNames = {};
config.useInequalityConstraints = true;

% Set model stiffness parameters
p_vec = getModelParameters(config, 0, 0);
sf = p_vec(1)*p_vec(2)*p_vec(3);
config.paramValues.k_h = 20/sf;
config.paramValues.k_k = 20/sf;

config.optConfig.tolerance = 2;
config.optConfig.print_level = 5;

%% Set up the parameter grid for continuation study
nGrid     = numel(config.cont.gridParamNames); % Total grid dimensions
gridMin   = config.cont.gridParamMin;  % Minimum grid values
gridMax   = config.cont.gridParamMax;  % Maximum grid values  
gridStep  = config.cont.gridParamStep; % Grid step sizes
gridNames = config.cont.gridParamNames; % Grid variable names

% Batch size for parallel processing:
% n_batch = config.cont.n_batch;
n_batch = 2;
% Every n_out iteration is visulized and saved
n_out   = config.cont.n_out;

%% Initialize continuation data
filenameMAT = config.cont.filenameMAT;        % Base filename for continuation data. If exist, continuation will continue from that file.
warmstartfileMAT = config.cont.warmStartMAT;  % Warm start data file
% Check if continuation data already exists
if isfile(['data\',filenameMAT, '.mat']) 
    extract_data(filenameMAT);
    traj = interpolate_Z2traj(config, z_final(1,:)');
else
    [z_initial, z_final, gridPts, cost, status, traj, n_z] = initializeFromWarmstart(warmstartfileMAT, gridMin, gridMax, gridStep, config);
    % Meaning of status:
    %  -1 = not processed but a solution at this grid point already exists and has a similar cost value (as measured by config.cont.costNew)
    %   0 = not processed
    %   1 = succesfully processed and best result
    %   2 = succesfully processed but better result existed (as measured by config.cont.costSlack)
    %   3 = failed to process, optimization returned error. 
end

%% Initialize dynamic plotting for monitoring progress
hPlotsStruct = init_dynamic_plots(gridPts, cost, status, gridNames, config);

%% Restart the parallel computing pool
% cd(fileparts(mfilename('fullpath')))
% delete(gcp('nocreate'));  % close any existing pool
% parpool('local');         % start a fresh pool
%% Main continuation loop - process solutions until all are handled
% As long as there are unprocessed rows, do:
lastSize = size(status,1); % Track size for save intervals
try
    while sum(status<=0)>0 % Continue while unprocessed solutions remain
        
        % Sort solutions by processing priority (cost + status weighting)
        [~,indices] = sort(cost + (status==-1)*max(cost) + (status>0)*10*max(cost));
    
        % Pre-allocate batch processing variables
        z_finalRun = zeros(min(sum(status<=0), n_batch), n_z);
        costRun = zeros(min(sum(status<=0), n_batch), 1);
        convergedRun = false(min(sum(status<=0), n_batch), 1);
        
        % Extract batch data for parallel processing
        z_parfor = z_initial(indices(1:min(sum(status<=0),n_batch)),:);
        grid_parfor = gridPts(indices(1:min(sum(status<=0),n_batch)), :);
        
        % Process batch in parallel
        % parfor
        parfor i = 1:min(sum(status<=0),n_batch)  % for debuggining you need to turn off parallel computing
            % Process them in the optimizer:
            [z_finalRun(i,:), costRun(i), return_status] = continuation_single_step(z_parfor(i,:)', grid_parfor(i,:), traj, config);
            convergedRun(i) = return_status=="Solve_Succeeded";
        end
        %%
        % Update continuation strategy with batch results
        [succStruct, z_final, status, cost] = continuation_strategy(gridPts, cost, status, indices, costRun, z_finalRun, convergedRun, z_final, n_batch, config);
        % Filter data to remove suboptimal and failed solutions
        idx_filter = ~(status > 1); % Keep status -1, 0, 1
        
        % Merge base data with new successors
        z_final   = [z_final(idx_filter,:); succStruct.z_final];
        z_initial = [z_initial(idx_filter,:); succStruct.z_init];
        status    = [status(idx_filter); succStruct.status];
        cost      = [cost(idx_filter); succStruct.cost];
        gridPts   = [gridPts(idx_filter,:); succStruct.gridPts];

        % Update progress plots
        update_dynamic_plots(gridPts, cost, status, config, hPlotsStruct);
    
        % Save intermediate results at specified intervals
        if size(status,1) - lastSize > n_out 
            lastSize = size(status,1);
            save_data([dated_filename, '_UNFINISHED'], 'o', ...
                      'status', status, 'cost', cost, 'gridPts', gridPts, 'config', config);
        end
    end
catch ME
    % Emergency handling for crashes
    crashTime = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n========== CONTINUATION CRASHED ==========\n');
    fprintf('Time of crash: %s\n', crashTime);
    fprintf('Error message: %s\n', ME.message);
    fprintf('==========================================\n');

    % Emergency save of current state
    fprintf('Emergency save after crash...\n');
    save_data([dated_filename, '_CRASH'], 'o', 'z_final', z_final, 'z_initial', z_initial, ...
              'status', status, 'cost', cost, 'gridPts', gridPts, 'config', config);
end

%% Final processing and saving of results
% Filter to keep only optimal solutions (status == 1)
z_final   = z_final(status==1, :);
z_initial = z_initial(status==1, :);
cost      = cost(status==1, :);
gridPts   = gridPts(status==1, :);
status    = status(status==1, :);

% Save final results
filename = save_data(dated_filename, 'v', 'z_final', z_final, 'z_initial', z_initial, ...
                     'status', status, 'cost', cost, 'gridPts', gridPts, 'config', config);
documentInequalityConstraints(config, [filename,'_ineq.txt']) % Document constraint definitions
 
%%
function hPlotsStruct = init_dynamic_plots(grid_array, cost, status, gridNames, config)
% Initialize dynamic plots for monitoring continuation progress

    views = config.cont.plotViews;        % Parameter pairs to visualize
    sliceVals = config.cont.sliceValues;  % Fixed values for higher dimensions
    nViews = numel(views);

    figure(1);
    clf;

    % Set up subplot grid
    nCols = ceil(sqrt(nViews));
    nRows = ceil(nViews / nCols);

    % Initialize plot handles structure
    hPlotsStruct = struct();
    hPlotsStruct.axes = gobjects(nViews,1);
    hPlotsStruct.hLater   = gobjects(nViews,1);
    hPlotsStruct.hWaiting = gobjects(nViews,1);
    hPlotsStruct.hValid   = gobjects(nViews,1);
    hPlotsStruct.hSubpar  = gobjects(nViews,1);
    hPlotsStruct.hFailed  = gobjects(nViews,1);

    % Create subplot for each view
    for i = 1:nViews
        dimsToPlot = views{i};  % Dimensions to plot on x and y axes
        tol = config.cont.sliceTolerance; % Tolerance for slice matching

        % Identify dimensions to keep fixed (slicing dimensions)
        dimsSlice = setdiff(1:size(grid_array,2), dimsToPlot);

        % Create mask for points matching slice values
        mask = true(size(grid_array,1),1);
        for d = dimsSlice
            mask = mask & abs(grid_array(:,d) - sliceVals(d)) < tol;
        end

        % Extract data for plotting
        x = grid_array(mask, dimsToPlot(1));
        y = grid_array(mask, dimsToPlot(2));
        z = cost(mask);
        status_vec = status(mask);

        % Create subplot with colored points by status
        hPlotsStruct.axes(i) = subplot(nRows, nCols, i);
        hold on; grid on; box on
        
        % Plot points by status with different colors/markers
        hPlotsStruct.hSubpar(i)  = safe_plot3(x(status_vec==2),  y(status_vec==2),  z(status_vec==2),  'k.');  % Suboptimal
        hPlotsStruct.hWaiting(i) = safe_plot3(x(status_vec==0),  y(status_vec==0),  z(status_vec==0),  'bo');  % Waiting
        hPlotsStruct.hValid(i)   = safe_plot3(x(status_vec==1),  y(status_vec==1),  z(status_vec==1),  'gs');  % Valid/Optimal
        hPlotsStruct.hLater(i)   = safe_plot3(x(status_vec==-1), y(status_vec==-1), z(status_vec==-1), 'm.'); % Postponed
        hPlotsStruct.hFailed(i)  = safe_plot3(x(status_vec==3),  y(status_vec==3),  z(status_vec==3),  'rx');  % Failed
        
        % Label axes and title
        xlabel(gridNames{dimsToPlot(1)});
        ylabel(gridNames{dimsToPlot(2)});
        zlabel('Cost');
        title(sprintf('View: Param %s vs Param %s', gridNames{dimsToPlot(1)}, gridNames{dimsToPlot(2)}));

        axis tight
        view(hPlotsStruct.axes(i), [-75 15]); % Set 3D view angle
        legend('Suboptimal','Waiting','Valid','Later','Failed','Location','bestoutside')
    end
end

%%
function update_dynamic_plots(grid_array, cost, status, config, hPlotsStruct)
% Update dynamic plots with current continuation state

    views = config.cont.plotViews;
    sliceVals = config.cont.sliceValues;
    tol = config.cont.sliceTolerance;
    nViews = numel(views);

    for i = 1:nViews
        dimsToPlot = views{i};
        
        % Identify slicing dimensions
        dimsSlice = setdiff(1:size(grid_array,2), dimsToPlot);
        
        % Create mask for current slice
        mask = true(size(grid_array,1),1);
        for d = dimsSlice
            mask = mask & abs(grid_array(:,d) - sliceVals(d)) < tol;
        end
        
        % Extract data for current slice
        x = grid_array(mask, dimsToPlot(1));
        y = grid_array(mask, dimsToPlot(2));
        z = cost(mask);
        status_vec = status(mask);

        % Update each status category plot
        safe_update_plot3(hPlotsStruct.hSubpar(i),  x, y, z, status_vec, 2);  % Suboptimal
        safe_update_plot3(hPlotsStruct.hWaiting(i), x, y, z, status_vec, 0);  % Waiting
        safe_update_plot3(hPlotsStruct.hValid(i),   x, y, z, status_vec, 1);  % Valid/Optimal
        safe_update_plot3(hPlotsStruct.hFailed(i),  x, y, z, status_vec, 3);  % Failed
        safe_update_plot3(hPlotsStruct.hLater(i),   x, y, z, status_vec, -1); % Postponed
    end

    drawnow limitrate % Update display efficiently
end

function h = safe_plot3(x, y, z, varargin)
% Safe plotting that handles empty data
    if isempty(x)
        h = plot3(nan, nan, nan, varargin{:}, 'Visible', 'off'); % Invisible placeholder
    else
        h = plot3(x, y, z, varargin{:}); % Normal plot
    end
end

function safe_update_plot3(h, x, y, z, status_vec, status_code)
% Safely update plot data based on status
    idx = (status_vec == status_code);
    if any(idx)
        set(h, 'XData', x(idx), 'YData', y(idx), 'ZData', z(idx), 'Visible', 'on');
    else
        set(h, 'XData', nan, 'YData', nan, 'ZData', nan, 'Visible', 'off'); % Hide if no data
    end
end

function [z_initial, z_final, gridPts, cost, status, traj, n_z] = initializeFromWarmstart(warmstartfileMAT, gridMin, gridMax, gridStep, config)
% INITIALIZEFROMWARMSTART Initialize continuation data from warm start results
%
%   Loads and filters warm start solutions, verifies them through parallel
%   optimization, and initializes continuation data structures.
%
%   Inputs:
%     warmstartfileMAT - filename of warm start data
%     gridMin, gridMax, gridStep - grid boundaries and step sizes
%     config - system configuration structure
%
%   Outputs:
%     z_initial - initial guess variables for continuation
%     z_final   - placeholder for final solution variables
%     gridPts   - grid points for continuation
%     cost      - cost values (initialized as infinite)
%     status    - status codes for each point (initialized as 0)
%     traj      - reference trajectory for interpolation

    % Verify warm start file exists
    if ~isfile(['data\', warmstartfileMAT, '.mat'])
        error('Warm start file does not exist: %s', warmstartfileMAT);
    end
    
    p_vec = getModelParameters(config, 0, 0);
    sf = p_vec(1)*p_vec(2)*p_vec(3);

    % Load warm start data
    load(warmstartfileMAT, 'Z_array', 'cost_array', 'gridPoint_array', 'numConverged_array', 'largestNorm');

    v_target    = config.operatingCond.v_avg;
    [~, optIdx] = min(abs(gridPoint_array(:,1) - v_target));

    n_z = size(Z_array,2);
    z_warmstart = Z_array(optIdx, :);
    traj = interpolate_Z2traj(config, Z_array(optIdx(1),:)');

    % Align grid points to grid
    gridPoint = gridPoint_array(optIdx,2:3);
    k = round(gridPoint./gridStep);
    gridPts = (k .* gridStep);
    
    % Verify warm start solutions through parallel optimization
    for i = 1:numel(optIdx)
        [z_finalRun(i,:), costRun(i), return_status] = continuation_single_step(z_warmstart(i,:)', gridPts(i,:), traj, config);
        convergedRun(i) = return_status=="Solve_Succeeded";
    end
    
    % Initialize continuation data structures
    cost = costRun;     % Costs initialized as infinite
    z_initial = z_finalRun; % Initial guesses from warm start
    z_final = zeros(numel(optIdx), n_z);      % Placeholder for final solutions
    status = zeros(numel(optIdx), 1);         % All points unprocessed (status 0)
end
