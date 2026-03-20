function plot_traj(config, input, dt, trajLB, trajUB)
% PLOT_TRAJ Plots state, input, and ground reaction force trajectories
%   Handles both full and reduced trajectories with optional bounds
%
%   Inputs:
%     config  - system configuration structure
%     input   - trajectory data (struct or optimization output)
%     dt      - time step for interpolation (optional)
%     trajLB  - lower bound trajectories (optional)
%     trajUB  - upper bound trajectories (optional)
%
%   author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    if nargin < 2
        error('Not enough input arguments');
    end
    
    % Process input data to full trajectory format
    if isstruct(input) && isfield(input, 'x')
        if size(input(1).x, 2) == length(config.defaultStateValues)
            traj_full = input;
            all_x = cat(1, traj_full(:).x);
            all_u = cat(1, traj_full(:).u);
            active_x = find(any(all_x ~= 0, 1));
            active_u = find(any(all_u ~= 0, 1));
        else
            traj_full = map2ClassTrajectory(config, input);
            all_x = cat(1, traj_full(:).x);
            all_u = cat(1, traj_full(:).u);
            get_indices = str2func([config.model_name, '.', config.phaseSequence{1}, '.get_classIndices']);
            [x_idx, u_idx] = get_indices();
            active_x = x_idx;
            active_u = u_idx;
        end
    else
        if nargin < 3
            dt = 0.01;
        end
        traj_interp = interpolate_Z2traj(config, input, [], 1, dt);
        traj_full = map2ClassTrajectory(config, traj_interp);
        all_x = cat(1, traj_full(:).x);
        all_u = cat(1, traj_full(:).u);
        get_indices = str2func([config.model_name, '.', config.phaseSequence{1}, '.get_classIndices']);
        [x_idx, u_idx] = get_indices();
        active_x = x_idx;
        active_u = u_idx;
    end

    % Process trajectory bounds if provided
    all_x_LB = []; all_u_LB = [];
    all_x_UB = []; all_u_UB = [];
    
    if nargin >= 4 && ~isempty(trajLB)
        trajLB_full = map2ClassTrajectory(config, trajLB);
        all_x_LB = cat(1, trajLB_full(:).x);
        all_u_LB = cat(1, trajLB_full(:).u);
    end
    if nargin >= 5 && ~isempty(trajUB)
        trajUB_full = map2ClassTrajectory(config, trajUB);
        all_x_UB = cat(1, trajUB_full(:).x);
        all_u_UB = cat(1, trajUB_full(:).u);
    end

    % Create global time vector across all phases
    phases = config.phaseSequence;
    nPhases = numel(phases);
    totalT = [];
    for i = 1:nPhases
        if ~isfield(traj_full(i), 't') || isempty(traj_full(i).t)
            error('Phase %d missing time vector', i);
        end
        if i == 1
            totalT = traj_full(i).t;
        else
            totalT = [totalT; traj_full(i).t + totalT(end)];
        end
    end
    
    % Generate phase colors and labels
    phase_colors = lines(nPhases);
    phase_labels = cell(nPhases, 1);
    for i = 1:nPhases
        phase_count = sum(strcmp(phases(1:i), phases{i}));
        phase_labels{i} = sprintf('%s %d', phases{i}, phase_count);
    end
    
    % Plot state trajectories
    fig_states = figure('Name', 'State Trajectories', 'NumberTitle', 'off');
    nX = numel(active_x);
    nCols = ceil(sqrt(nX));
    nRows = ceil(nX / nCols);
    
    legend_handles = cell(nPhases, 1);
    for k = 1:nX
        subplot(nRows, nCols, k);
        hold on;
        idxStart = 1;
        for i = 1:nPhases
            N = size(traj_full(i).x, 1);
            idxEnd = idxStart + N - 1;
            get_indices = str2func([config.model_name, '.', config.phaseSequence{i}, '.get_classIndices']);
            [IndicesPerPhase_x, ~] = get_indices();
            
            % Plot state trajectory
            h = plot(totalT(idxStart:idxEnd), all_x(idxStart:idxEnd, active_x(k)), ...
                    'Color', phase_colors(i,:), 'LineWidth', 1.5);
            if k == 1
                legend_handles{i} = h;
            end
            
            % Plot bounds if available
            plotBounds(totalT(idxStart:idxEnd), all_x_LB, all_x_UB, ...
                      idxStart, idxEnd, active_x(k), IndicesPerPhase_x);
            
            idxStart = idxEnd + 1;
        end
        title(config.defaultStateNames{active_x(k)}, 'Interpreter', 'none');
        xlabel('Time [s]');
        grid on;
    end
    legend([legend_handles{:}], phase_labels, 'Location', 'bestoutside');
    
    % Plot input trajectories
    fig_inputs = figure('Name', 'Input Trajectories', 'NumberTitle', 'off');
    nU = numel(active_u);
    nCols = ceil(sqrt(nU));
    nRows = ceil(nU / nCols);
    
    for k = 1:nU
        subplot(nRows, nCols, k);
        hold on;
        idxStart = 1;
        for i = 1:nPhases
            N = size(traj_full(i).u, 1);
            idxEnd = idxStart + N - 1;
            get_indices = str2func([config.model_name, '.', config.phaseSequence{i}, '.get_classIndices']);
            [~, IndicesPerPhase_u] = get_indices();
            
            plot(totalT(idxStart:idxEnd), all_u(idxStart:idxEnd, active_u(k)), ...
                 'Color', phase_colors(i,:), 'LineWidth', 1.5);

            % Plot bounds if available
            plotBounds(totalT(idxStart:idxEnd), all_u_LB, all_u_UB, ...
                      idxStart, idxEnd, active_u(k), IndicesPerPhase_u);
            
            idxStart = idxEnd + 1;
        end
        title(config.defaultInputNames{active_u(k)}, 'Interpreter', 'none');
        xlabel('Time [s]');
        grid on;
    end
    legend(phase_labels, 'Location', 'bestoutside');
    
    % Plot ground reaction forces
    plotGroundReactionForces(traj_full, phases, phase_colors, phase_labels, totalT);
end

function plotBounds(time_segment, data_LB, data_UB, idxStart, idxEnd, active_idx, valid_indices)
% PLOTBOUNDS Helper function to plot trajectory bounds
    if ~isempty(data_LB) && ~isempty(data_UB) && any(active_idx == valid_indices)
        lb_vals = data_LB(idxStart:idxEnd, active_idx);
        ub_vals = data_UB(idxStart:idxEnd, active_idx);
        
        if any(abs(lb_vals - ub_vals) > 1e-5)
            plot(time_segment, lb_vals, 'k--', 'LineWidth', 1);
            plot(time_segment, ub_vals, 'k--', 'LineWidth', 1);
        end
    else
        if ~isempty(data_LB) && any(active_idx == valid_indices)
            plot(time_segment, data_LB(idxStart:idxEnd, active_idx), 'k--', 'LineWidth', 1);
        end
        if ~isempty(data_UB) && any(active_idx == valid_indices)
            plot(time_segment, data_UB(idxStart:idxEnd, active_idx), 'k--', 'LineWidth', 1);
        end
    end
end

function plotGroundReactionForces(traj_full, phases, phase_colors, phase_labels, totalT)
% PLOTGROUNDREACTIONFORCES Plot ground reaction force components
    nPhases = numel(phases);
    
    % Determine maximum number of lambda components across all phases
    max_lambda_components = 0;
    for i = 1:nPhases
        current_lambda = traj_full(i).lambda;
        if isrow(current_lambda)
            current_lambda = current_lambda';
        end
        max_lambda_components = max(max_lambda_components, size(current_lambda, 2));
    end
    
    % Generate component names based on number of components
    if max_lambda_components == 1
        comp_names = {'Vertical GRF'};
    elseif max_lambda_components == 2
        comp_names = {'Horizontal GRF', 'Vertical GRF'};
    elseif max_lambda_components == 3
        comp_names = {'X GRF', 'Y GRF', 'Z GRF'};
    else
        comp_names = arrayfun(@(x) sprintf('GRF Component %d', x), 1:max_lambda_components, 'UniformOutput', false);
    end
    
    % Create figure and determine subplot layout
    fig_lambda = figure('Name', 'Ground Reaction Forces', 'NumberTitle', 'off');
    
    if max_lambda_components > 3
        nCols = ceil(sqrt(max_lambda_components));
        nRows = ceil(max_lambda_components / nCols);
    else
        nCols = max_lambda_components;
        nRows = 1;
    end
    
    % Plot each lambda component
    for k = 1:max_lambda_components
        if max_lambda_components > 1
            subplot(nRows, nCols, k);
        end
        hold on;
        idxStart = 1;
        for i = 1:nPhases
            N = size(traj_full(i).x, 1);
            idxEnd = idxStart + N - 1;
            
            current_lambda = traj_full(i).lambda;
            if isrow(current_lambda)
                current_lambda = current_lambda';
            end
            
            % Plot component if it exists for this phase
            if size(current_lambda, 2) >= k
                plot(totalT(idxStart:idxEnd), current_lambda(1:N, k), ...
                     'Color', phase_colors(i,:), 'LineWidth', 1.5);
            end
            idxStart = idxEnd + 1;
        end
        title(comp_names{k});
        xlabel('Time [s]');
        grid on;
    end
    legend(phase_labels, 'Location', 'bestoutside');
end