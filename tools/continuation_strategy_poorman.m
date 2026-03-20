function [succStruct, z_final, status, cost] = continuation_strategy_poorman(grid_array, cost, status, indices, costRun, z_finalRun, convergedRun, z_final, n_batch, config)
% CONTINUATION_STRATEGY Manage grid-based continuation strategy for parameter studies
%
%   Implements a systematic approach for exploring parameter spaces by
%   propagating successful solutions to neighboring grid points. Handles
%   solution ranking, successor generation, and grid boundary management.
%
%   Inputs:
%     grid_array   - M×N matrix of grid points (M points, N parameters)
%     cost         - M×1 vector of cost values for each grid point
%     status       - M×1 vector of status codes for each grid point
%     indices      - indices of current batch being processed
%     costRun      - cost values from current optimization batch
%     z_finalRun   - decision variables from current optimization batch
%     convergedRun - convergence flags from current optimization batch
%     z_final      - M×n_z matrix of final decision variables
%     n_batch      - number of points in current batch
%     config       - configuration structure with continuation parameters
%
%   Outputs:
%     succStruct   - structure containing successor information:
%                    status: status codes for successors (-1=similar exists, 0=new)
%                    gridPts: grid coordinates of successors
%                    cost: cost values from parent solutions
%                    z_init: initial guesses for successors (parent solutions)
%                    z_final: placeholder for successor solutions
%     z_final      - updated decision variable matrix
%     status       - updated status vector
%     cost         - updated cost vector
%
%   Status Codes:
%     0: Untried   1: Optimal   2: Suboptimal   3: Failed
%     -1: Successor with similar solution exists
%
% author: C. David Remy, Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025


    n_gridDim = size(grid_array, 2);    % Number of grid parameters (e.g., 2 for v,k)
    n_z = size(z_finalRun, 2);          % Length of z vectors
    
    % Determine number of successful solutions in current batch
    n_success = min(sum(status(indices) <= 0), n_batch);
    succ_alloc = 2 * n_gridDim * n_success; % Allocate for 2 successors per dimension
    
    % Preallocate successor structure
    succStruct.status    = nan(succ_alloc,1);
    succStruct.gridPts   = nan(succ_alloc, n_gridDim);
    succStruct.cost      = nan(succ_alloc,1);
    succStruct.z_init    = nan(succ_alloc, n_z);
    succStruct.z_final   = nan(succ_alloc, n_z);
    
    succCounter = 0;

    % Extract grid boundaries and step sizes from configuration.
    % Use integer grid indices to avoid float-equality issues.
    grid_min  = config.cont.gridParamMin;
    grid_max  = config.cont.gridParamMax;
    grid_step = config.cont.gridParamStep;
    grid_min  = grid_min(:).';
    grid_max  = grid_max(:).';
    grid_step = grid_step(:).';
    grid_idx_all = round((grid_array - grid_min) ./ grid_step);
    grid_idx_min = zeros(1, n_gridDim);
    grid_idx_max = round((grid_max - grid_min) ./ grid_step);
        
    % Process each solution in the current batch
    for i = 1:n_success
        idx = indices(i);
    
        % Update cost and solution for current grid point
        cost(idx) = costRun(i);
        z_final(idx,:) = z_finalRun(i,:);
    
        if ~convergedRun(i)
            status(idx) = 3; % Mark as failed
            continue;
        end
    
        % Check if better solution already exists at this grid point
        grid_idx = grid_idx_all(idx, :);
        same_point = all(grid_idx_all == grid_idx, 2);
        betterExists = any(cost(same_point & status==1) < cost(idx) + config.cont.costSlack);
    
        if betterExists
            status(idx) = 2; % Mark as suboptimal (worse than existing)
            continue;
        end
    
        % This solution is the new best at this grid point
        status(idx) = 1; % Mark as optimal
        same_point(idx) = false; % Exclude current index from same_point
        status(same_point & status==1) = 2; % Demote previous optimal solutions
    
        % Generate successors in each grid dimension
        for d = 1:n_gridDim
            for direction = [-1, 1] % Explore both directions along each dimension
                % Calculate neighbor grid index
                new_idx = grid_idx(d) + direction;

                % Check if new index is within grid boundaries
                if new_idx < grid_idx_min(d) || new_idx > grid_idx_max(d)
                    continue; % Skip if outside grid
                end
    
                % Construct target grid point
                grid_target_idx = grid_idx;
                grid_target_idx(d) = new_idx;
                grid_target = grid_min + grid_target_idx .* grid_step;
    
                % Check if similar solution already exists at target
                match = all(grid_idx_all == grid_target_idx, 2) & status ~= 3;
                similarExists = any(cost(match) < cost(idx) + config.cont.costNew);
    
                % Add successor to structure
                succCounter = succCounter + 1;
                succStruct.status(succCounter)    = ternary(similarExists, -1, 0);
                succStruct.gridPts(succCounter,:) = grid_target;
                succStruct.cost(succCounter)      = cost(idx);
                succStruct.z_init(succCounter,:)  = z_final(idx,:);
                succStruct.z_final(succCounter,:) = zeros(1, n_z); % Placeholder
            end
        end
    end
    
    % Trim unused preallocated entries
    succStruct.status    = succStruct.status(1:succCounter);
    succStruct.gridPts   = succStruct.gridPts(1:succCounter,:);
    succStruct.cost      = succStruct.cost(1:succCounter);
    succStruct.z_init    = succStruct.z_init(1:succCounter,:);
    succStruct.z_final   = succStruct.z_final(1:succCounter,:);
end

function val = ternary(cond, valTrue, valFalse)
% TERNARY Ternary conditional operator
%
%   Provides ternary operator functionality: returns valTrue if cond is true,
%   otherwise returns valFalse.
%
%   Inputs:
%     cond    - boolean condition
%     valTrue - value to return if condition is true
%     valFalse- value to return if condition is false
%
%   Output:
%     val     - either valTrue or valFalse based on condition

    if cond
        val = valTrue;
    else
        val = valFalse;
    end
end
