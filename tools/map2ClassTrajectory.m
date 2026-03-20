function traj_full = map2ClassTrajectory(config, traj)
% MAP2CLASSTRAJECTORY Maps reduced trajectories to full states using relations
%   Handles phase transitions and maintains kinematic consistency
%   Inputs:
%     config - system configuration structure
%     traj   - trajectory structure containing states, controls, time, 
%              parameters, and constraint forces in the reduced coordiante
%              form
%
%   Output:
%     traj   - trajectory structure containing the FULL states, controls, time, 
%              parameters, and constraint of the specific class. 
% The class is specified by the fisrt letter of the config.model_name,
% e.g., HF_min --> Hopper
%    
%    author: Iskandar Khemakhem, Oussama Barhoumi, IAMS, Uni Stuttgart, 2025

    % Get parameter functions
    getParameterNames = str2func([config.model_name, '.getParameterNames']);
    
    % Verify model and initialize
    if startsWith(config.model_name, 'h', 'IgnoreCase', true)
        % Hopper model parameters
        modelTarget = 'hopper5DOF';
        nx_full = 10;
        nu_full = 2; 
        x_temp = zeros(nx_full, 1);

        nPhases = numel(traj);
        traj_full = traj;  % Preallocate structure
        paramNames = getParameterNames(config);
        params = getModelParameters(config, 0, 0);
        
        % Phi angle handling
        if ismember('phi_0', paramNames)
            phi_fun = @(params, config) params(getParameterIndex('phi_0', config));
        else
            phi_fun = @(~, ~) 0;
        end

        % Kinematic relations for missing states
        relations = {
            {'1', @(x, params, config, x_temp) ...
                x_temp(1) + ...
                (phi_fun(params, config) + x_temp(4)) * params(getParameterIndex('r_f', config)) + ...
                x_temp(5) * sin(phi_fun(params, config) + x_temp(4)) - ...
                params(getParameterIndex('r_f', config)) * (phi_fun(params, config) + x(:,4)) - ...
                x(:,5) .* sin(phi_fun(params, config) + x(:,4))};
        
            {'2', @(x, params, config, x_temp) ...
                x(:,5) .* cos(x(:,4) + phi_fun(params, config)) + ...
                params(getParameterIndex('r_f', config))};
        
            {'6', @(x, params, config, x_temp) ...
                -params(getParameterIndex('r_f', config)) * x(:,9) - ...
                x(:,10) .* sin(phi_fun(params, config) + x(:,4)) - ...
                x(:,5) .* x(:,9) .* cos(phi_fun(params, config) + x(:,4))};
        
            {'7', @(x, params, config, x_temp) ...
                x(:,10) .* cos(phi_fun(params, config) + x(:,4)) - ...
                x(:,5) .* x(:,9) .* sin(phi_fun(params, config) + x(:,4))}
        };
        
    elseif startsWith(config.model_name, 'b', 'IgnoreCase', true)
        % Hopper model parameters
        modelTarget = 'biped5DOF';
        nx_full = 10;
        nu_full = 4; 
        x_temp = zeros(nx_full, 1);

        nPhases = numel(traj);
        traj_full = traj;  % Preallocate structure
        paramNames = getParameterNames(config);
        params = getModelParameters(config, 0, 0);
        
        % Phi angle     
        if ismember('phi_0', paramNames)
            phi_fun = @(params, config) params(getParameterIndex('phi_0', config));
        else
            phi_fun = @(~, ~) 0;
        end

        % Kinematic relations for missing states
        relations = {
            {'1', @(x, params, config, x_temp) ...
                x_temp(1) + ...
                (phi_fun(params, config) + x_temp(4)) * params(getParameterIndex('r_f', config)) + ...
                x_temp(5) * sin(phi_fun(params, config) + x_temp(4)) - ...
                params(getParameterIndex('r_f', config)) * (phi_fun(params, config) + x(:,4)) - ...
                x(:,5) .* sin(phi_fun(params, config) + x(:,4))};
        
            {'2', @(x, params, config, x_temp) ...
                x(:,5) .* cos(x(:,4) + phi_fun(params, config)) + ...
                params(getParameterIndex('r_f', config))};
        
            {'6', @(x, params, config, x_temp) ...
                -params(getParameterIndex('r_f', config)) * x(:,9) - ...
                x(:,10) .* sin(phi_fun(params, config) + x(:,4)) - ...
                x(:,5) .* x(:,9) .* cos(phi_fun(params, config) + x(:,4))};
        
            {'7', @(x, params, config, x_temp) ...
                x(:,10) .* cos(phi_fun(params, config) + x(:,4)) - ...
                x(:,5) .* x(:,9) .* sin(phi_fun(params, config) + x(:,4))}
        };
    else
        error('Class not supported yet. Please adapt code to accomodate your class');
    end

    % Track previous phase's full state for continuity
    prev_full_x = config.defaultStateValues(:)';
    
    for i = 1:nPhases
        phase_name = config.phaseSequence{i};
        
        % Get optimized indices for current phase
        get_indices = str2func([config.model_name, '.', phase_name, '.get_classIndices']);
        [x_idx, u_idx] = get_indices();

        % Verify dimensions
        if numel(config.defaultStateValues) ~= nx_full
            error('defaultStateValues must be of length %d', nx_full);
        end
        if numel(config.defaultInputValues) ~= nu_full
            error('defaultInputValues must be of length %d', nu_full);
        end

        % number of segments
        N = size(traj(i).x, 1);
        
        % Initialize with previous phase's states (for continuity)
        traj_full(i).x = repmat(prev_full_x, N, 1);

        % when entering a Stance for _min models, lock reference from end of previous Flight
        if strcmp(config.model_name(end-3:end), '_min') 
            if startsWith(config.model_name, 'h', 'IgnoreCase', true) && strcmp(phase_name, 'Stance')
                if i == 1
                    error('HF_min cannot start with Stance: there is no preceding Flight to lock values from.');
                end
                % Prefer explicit check that previous was Flight; if not, still use prev_full_x but warn
                if ~strcmp(config.phaseSequence{i-1}, 'Flight')
                    warning('Entering Stance but previous phase is not Flight; using last full state as reference.');
                end
                % prev_full_x should contain the end full-state of the previous phase (Flight)
                x_temp(1) = prev_full_x(1);
                x_temp(4) = prev_full_x(4);
                x_temp(5) = prev_full_x(5);
            elseif ~startsWith(config.model_name, 'h', 'IgnoreCase', true)
                warning('treating minimal models for your class is missing')
            end
        end
        
        % Fill in the optimized states with correct indexing
        for k = 1:length(x_idx)
            traj_full(i).x(:, x_idx(k)) = traj(i).x(:, k);
        end
        
        % Apply kinematic relations for missing states
        for j = 1:length(relations)
            relation = relations{j};
            stateIdx = str2double(relation{1});
            relationFunc = relation{2};
            
            if ~ismember(stateIdx, x_idx)
                traj_full(i).x(:, stateIdx) = relationFunc(traj_full(i).x, params, config, x_temp);
            end
        end
        
        % Handle inputs
        traj_full(i).u = repmat(config.defaultInputValues(:)', N, 1);
        for k = 1:length(u_idx)
            traj_full(i).u(:, u_idx(k)) = traj(i).u(:, k);
        end
        
        % Store for next phase initialization
        prev_full_x = traj_full(i).x(end, :);
        
        % Copy time vector
        if isfield(traj(i), 't')
            traj_full(i).t = traj(i).t;
        end
    end
end