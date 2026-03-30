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
dated_filename = append_date('PC_kh');
% % main continutation file
% configfile_HF_noPitch_min_Mac;
% % save initial config
% config_0 = config;
% % change some settings in the config
% config.optParameterNames = {};
% config.optParameterGlobal = [];
% config.optParameterConstant = [];
% config.optParameterInit = [];
% config.optParameterLowerBound = [];
% config.optParameterUpperBound = [];
% config.useInequalityConstraints = true;
% 
% %% set fixed value
% % config.paramValues.k_h = 1;
% config.paramValues.k_l = 200;

%% Restart the parallel computing pool
cd(fileparts(mfilename('fullpath')))
delete(gcp('nocreate'));  % close any existing pool
parpool('local');         % start a fresh pool
%% Set up the parameter grid for continuation study
nGrid     = numel(config.cont.gridOpCondNames) + numel(config.cont.gridParamNames); % Total grid dimensions
gridMin   = [config.cont.gridOpCondMin; config.cont.gridParamMin];  % Minimum grid values
gridMax   = [config.cont.gridOpCondMax; config.cont.gridParamMax];  % Maximum grid values  
gridStep  = [config.cont.gridOpCondStep; config.cont.gridParamStep]; % Grid step sizes
gridNames = {config.cont.gridOpCondNames; config.cont.gridParamNames}; % Grid variable names

% Batch size for parallel processing:
n_batch = config.cont.n_batch;
% Every n_out iteration is visulized and saved
n_out   = config.cont.n_out;

%% Initialize continuation data
filenameMAT = config.cont.filenameMAT;        % Base filename for continuation data. If exist, continuation will continue from that file.
warmstartfileMAT = config.cont.warmStartMAT;  % Warm start data file
% Check if continuation data already exists
if isfile([filenameMAT, '.mat']) 
    extract_data(filenameMAT);
    traj = interpolate_Z2traj(config, z_final(1,:)');
    if exist('z_initial', 'var') && ~isempty(z_initial)
        n_z = size(z_initial, 2);
    else
        n_z = size(z_final, 2);
    end
else
    % Original warm-start handling:
    % [z_initial, z_final, gridPts, cost, status, traj, n_z] = initializeFromWarmstart(warmstartfileMAT, gridMin, gridMax, gridStep, config);
    % Robust warm-start handling (works with warm-start maps from one or multiple optimized parameters):
    [z_initial, z_final, gridPts, cost, status, traj, n_z] = initializeFromWarmstart_robust(warmstartfileMAT, gridMin, gridMax, gridStep, config);
    % Meaning ofstatus: 
    %  -1 = not processed but a solution at this grid point already exists and has a similar cost value (as measured by config.cont.costNew)
    %   0 = not processed
    %   1 = succesfully processed and best result
    %   2 = succesfully processed but better result existed (as measured by config.cont.costSlack)
    %   3 = failed to process, optimization returned error. 
end

%% Initialize dynamic plotting for monitoring progress
hPlotsStruct = init_dynamic_plots(gridPts, cost, status, gridNames, config);

%% Main continuation loop - process solutions until all are handled
% As long as there are unprocessed rows, do:
lastSize = size(status,1); % Track size for save intervals
try
    while sum(status<=0)>0 % Continue while unprocessed solutions remain
        
        % Sort solutions by processing priority (finite-safe weighting).
        finite_cost = cost(isfinite(cost));
        if isempty(finite_cost)
            finite_max = 1;
        else
            finite_max = max(finite_cost);
            if finite_max <= 0
                finite_max = 1;
            end
        end

        priority = cost;
        priority(~isfinite(priority)) = finite_max + 1;
        priority(status==-1) = priority(status==-1) + (finite_max + 1);
        priority(status>0)   = priority(status>0)   + 2*(finite_max + 1);
        [~,indices] = sort(priority, 'ascend');
        n_proc = min(sum(status<=0), n_batch);
    
        % Pre-allocate batch processing variables
        z_finalRun   = zeros(n_proc, n_z);
        costRun      = zeros(n_proc, 1);
        convergedRun = false(n_proc, 1);
        
        % Extract batch data for parallel processing
        z_parfor = z_initial(indices(1:n_proc),:);
        grid_parfor = gridPts(indices(1:n_proc), :);
        
        % Process batch in parallel
        parfor i = 1:n_proc
        %for i = 1:n_proc  % for debugging you can keep this as a regular for-loop
            % Process them in the optimizer:
            [z_finalRun(i,:), costRun(i), return_status] = continuation_single_step(z_parfor(i,:)', grid_parfor(i,:), traj, config);
            convergedRun(i) = return_status=="Solve_Succeeded";
        end
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
        xlabel(gridNames{dimsToPlot(1)}{1});
        ylabel(gridNames{dimsToPlot(2)}{1});
        zlabel('Cost');
        title(sprintf('View: Param %s vs Param %s', gridNames{dimsToPlot(1)}{1}, gridNames{dimsToPlot(2)}{1}));

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

    warmstartPath = fullfile('CoCoLeg_updated', 'data', [warmstartfileMAT, '.mat']);

    % Verify warm start file exists
    if ~isfile(warmstartPath)
        error('Warm start file does not exist: %s', warmstartfileMAT);
    end
    
    % Load warm start data
    load(warmstartPath, 'Z_array', 'cost_array', 'gridPoint_array', 'numConverged_array', 'largestNorm');
    
    % Filter solutions using moving median outlier detection
    window = 3; % Number of neighbors for median filter
    localMed = movmedian(cost_array, 2*window+1);
    jumpThreshold = 0.02; % Cost jump threshold for outlier detection
    isOutlier = (cost_array - localMed > jumpThreshold) & (cost_array ~= inf);
    isOutlier(1:3) = false; % Exclude first 3 points
    
    % Select optimal points meeting convergence and boundary criteria
    optIdx = find(cost_array ~= inf & numConverged_array >= 15 & ...
                  gridPoint_array(:,1) <= gridMax(1) & gridPoint_array(:,2) <= gridMax(2) & ...
                  gridPoint_array(:,1) >= gridMin(1) & gridPoint_array(:,2) >= gridMin(2) & ~isOutlier);

    if isempty(optIdx)
        error('No valid warm-start seeds found inside requested bounds.');
    end

    traj = interpolate_Z2traj(config, Z_array(optIdx(1),:)');
    z_ref = interpolate_traj2Z(config, traj);
    n_z = numel(z_ref);

    z_warmstart = zeros(numel(optIdx), n_z);
    z_warmstart(1,:) = z_ref(:).';
    for i = 2:numel(optIdx)
        traj_i = interpolate_Z2traj(config, Z_array(optIdx(i),:)');
        z_warmstart(i,:) = interpolate_traj2Z(config, traj_i).';
    end
    
    % Align grid points to grid
    k = round(gridPoint_array(optIdx,:)'./gridStep);
    gridPts = (k .* gridStep)';
    
    % Verify warm start solutions through parallel optimization
    z_finalRun = zeros(numel(optIdx), n_z);
    costRun = inf(numel(optIdx), 1);
    convergedRun = false(numel(optIdx), 1); 
    parfor i = 1:numel(optIdx)
        [z_finalRun(i,:), costRun(i), return_status] = continuation_single_step(z_warmstart(i,:)', gridPts(i,:), traj, config);
        convergedRun(i) = return_status=="Solve_Succeeded";
    end
    
    % Filter out solutions that degraded significantly
    badSeedsIdx = abs(costRun(:) - cost_array(optIdx)) > 0.02;
    optIdx = optIdx(~badSeedsIdx);
    gridPts = gridPts(~badSeedsIdx, :);
    
    % Initialize continuation data structures
    cost = inf(numel(optIdx), 1);     % Costs initialized as infinite
    z_initial = z_warmstart(~badSeedsIdx, :); % Initial guesses from warm start
    z_final = zeros(numel(optIdx), n_z);      % Placeholder for final solutions
    status = zeros(numel(optIdx), 1);         % All points unprocessed (status 0)
end

function [z_initial, z_final, gridPts, cost, status, traj, n_z] = initializeFromWarmstart_robust(warmstartfileMAT, gridMin, gridMax, gridStep, config)
% INITIALIZEFROMWARMSTART_ROBUST
% Robust warm-start initialization for continuation.
% Supports warm-start maps generated with one optimized parameter or multiple
% optimized parameters. If warm-start includes extra optimized parameters that
% are fixed in continuation, seeds are filtered to stay close to the fixed
% values in config.paramValues.

    warmstartPath = fullfile('CoCoLeg_updated', 'data', [warmstartfileMAT, '.mat']);
    if ~isfile(warmstartPath)
        error('Warm start file does not exist: %s', warmstartfileMAT);
    end

    ws = load(warmstartPath, 'Z_array', 'cost_array', 'gridPoint_array', 'numConverged_array', 'largestNorm', 'config');
    Z_array = ws.Z_array;
    cost_array = ws.cost_array;
    gridPoint_array = ws.gridPoint_array;
    numConverged_array = ws.numConverged_array;

    if isfield(ws, 'config')
        warmConfig = ws.config;
    else
        warmConfig = struct();
    end

    % Outlier filtering
    window = 3;
    localMed = movmedian(cost_array, 2*window+1);
    jumpThreshold = 0.02;
    isOutlier = (cost_array - localMed > jumpThreshold) & (cost_array ~= inf);
    isOutlier(1:3) = false;

    % Continuation dimensions (e.g., [v_avg, k_h] or [v_avg, k_l]).
    contOpNames    = config.cont.gridOpCondNames(:);
    contParamNames = config.cont.gridParamNames(:);
    contNames      = [contOpNames; contParamNames];

    % Warm-start grid columns are [v_avg, optimized parameters...].
    if isfield(warmConfig, 'optParameterNames') && ~isempty(warmConfig.optParameterNames)
        warmOptNames = warmConfig.optParameterNames(:);
    else
        warmOptNames = {};
    end

    if isempty(contOpNames)
        opName = 'v_avg';
    else
        opName = contOpNames{1};
    end
    warmColNames = [{opName}; warmOptNames];

    % If warm-start metadata is missing, allow direct positional match.
    if numel(warmColNames) ~= size(gridPoint_array, 2)
        if size(gridPoint_array, 2) == numel(contNames)
            warmColNames = contNames;
        else
            error('Warm-start gridPoint_array has %d columns, cannot map to continuation dimensions.', ...
                  size(gridPoint_array, 2));
        end
    end

    % Map warm-start grid columns to continuation dimensions.
    gridPoint_cont = nan(size(gridPoint_array, 1), numel(contNames));
    for d = 1:numel(contNames)
        idxCol = find(strcmp(contNames{d}, warmColNames), 1);
        if isempty(idxCol)
            error('Continuation dimension "%s" not found in warm-start grid columns.', contNames{d});
        end
        gridPoint_cont(:, d) = gridPoint_array(:, idxCol);
    end

    % Keep only seeds close to fixed parameter values for parameters that are
    % present in warm-start but not part of continuation grid parameters.
    fixedMask = true(size(gridPoint_array, 1), 1);
    fixedParamNames = setdiff(warmOptNames, contParamNames, 'stable');
    for p = 1:numel(fixedParamNames)
        pName = fixedParamNames{p};
        idxCol = find(strcmp(pName, warmColNames), 1);
        if isempty(idxCol) || ~isfield(config, 'paramValues') || ~isfield(config.paramValues, pName)
            continue;
        end

        pTarget = config.paramValues.(pName);
        pTol = max(1e-6, 0.1 * max(abs(pTarget), 1));  % 10% tolerance
        fixedMask = fixedMask & abs(gridPoint_array(:, idxCol) - pTarget) <= pTol;
    end

    withinGrid = all(gridPoint_cont >= (gridMin(:).' - 10*eps) & ...
                     gridPoint_cont <= (gridMax(:).' + 10*eps), 2);

    optIdx = find(cost_array ~= inf & numConverged_array >= 10 & ...
                  withinGrid & fixedMask & ~isOutlier);

    if isempty(optIdx)
        error('No valid warm-start seeds found inside requested bounds (after fixed-parameter filtering).');
    end

    % Adapt warm-start decision vectors to current continuation setup.
    z_warmstart = [];
    n_z = 0;
    for i = 1:numel(optIdx)
        cfg_i = config;

        % Use warm-start optimized parameter values when reconstructing the
        % trajectory, so slack/global-dependent quantities stay consistent.
        for p = 1:numel(warmOptNames)
            pName = warmOptNames{p};
            idxCol = find(strcmp(pName, warmColNames), 1);
            if ~isempty(idxCol) 
                cfg_i.paramValues.(pName) = gridPoint_array(optIdx(i), idxCol);
            end
        end

        traj_i = interpolate_Z2traj(cfg_i, Z_array(optIdx(i),:)');
        z_i    = interpolate_traj2Z(cfg_i, traj_i);
        z_i    = z_i(:).';

        if i == 1
            n_z = numel(z_i);
            z_warmstart = zeros(numel(optIdx), n_z);
        end
        z_warmstart(i, :) = z_i;
    end

    traj = interpolate_Z2traj(config, z_warmstart(1,:)');

    % Build continuation grid points only in continuation dimensions.
    k = round((gridPoint_cont(optIdx, :) - gridMin(:).') ./ gridStep(:).');
    gridPts = gridMin(:).' + k .* gridStep(:).';

    % Verify warm-start seeds through one continuation step.
    z_finalRun = zeros(numel(optIdx), n_z);
    costRun    = inf(numel(optIdx), 1);
    convergedRun = false(numel(optIdx), 1);
    parfor i = 1:numel(optIdx)
        [z_finalRun(i,:), costRun(i), return_status] = continuation_single_step(z_warmstart(i,:)', gridPts(i,:), traj, config);
        convergedRun(i) = return_status=="Solve_Succeeded";
    end

    badSeedsIdx = abs(costRun(:) - cost_array(optIdx)) > 0.02;
    optIdx = optIdx(~badSeedsIdx);
    gridPts = gridPts(~badSeedsIdx, :);

    cost = inf(numel(optIdx), 1);
    z_initial = z_warmstart(~badSeedsIdx, :);
    z_final = zeros(numel(optIdx), n_z);
    status = zeros(numel(optIdx), 1);
end

