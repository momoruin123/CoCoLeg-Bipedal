function [activeBounds, trajLB, trajUB] = checkActiveIneqConstraints(traj, config, ineqTxtFile)
% CHECKACTIVEINEQCONSTRAINTS Identify active inequality constraints in trajectory
%
%   Analyzes a trajectory to determine which inequality constraints are active
%   (touching their bounds) and extracts the corresponding bound values for
%   visualization or analysis purposes.
%
%   Inputs:
%     traj        - trajectory structure array (one per phase)
%     config      - system configuration structure
%     ineqTxtFile - optional text file containing custom inequality constraints
%
%   Outputs:
%     activeBounds - structure array of active constraints with phase, name, and indices
%     trajLB      - lower bound values for each phase and variable type
%     trajUB      - upper bound values for each phase and variable type
%
%   author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    % Extract model information and initialize parameters
    modelName     = config.model_name;
    phaseSequence = config.phaseSequence;
    P_sym = getModelParameters(config, 1, 1);  % Symbolic parameters
    P_num = getModelParameters(config, 0, 0);  % Numerical parameters
    tol           = 1e-4; % Tolerance for considering a constraint active
    
    % Initialize output structures
    activeBounds = struct('phase', {}, 'name', {}, 'idx', {});
    trajLB = struct([]);
    trajUB = struct([]);
    
    % Load ground reaction force function
    grf_func = str2func([modelName, config.GRFfunc]); 
    
    % === Process external inequality constraint definitions if provided ===
    useExtr = (nargin >= 3 && ~isempty(ineqTxtFile) && ...
               isfile(['C:\Users\Iskandar\Desktop\Research\LegStiffness\git\data\',ineqTxtFile,'.txt']));
    
    if useExtr
        extrDir = fullfile('C:\Users\Iskandar\Desktop\Research\LegStiffness\git\data\', 'extr');
        txtFile = ['C:\Users\Iskandar\Desktop\Research\LegStiffness\git\data\',ineqTxtFile,'.txt'];
        
        % Create extraction directory if it doesn't exist
        if ~exist(extrDir, 'dir'), mkdir(extrDir); end
        
        % Read and parse the constraint definition file
        txt = fileread(txtFile);
        parts = regexp(txt, '########## Phase:\s*(\w+)\s*##########', 'split');
        names = regexp(txt, '########## Phase:\s*(\w+)\s*##########', 'tokens');
        
        % Create temporary model name to avoid namespace conflicts
        origModelName = modelName; % Preserve original model name
        modelName = [modelName, '_extr'];
        
        % Process each phase's constraint definitions
        for i = 1:length(names)
            phaseName = strtrim(names{i}{1});
            phaseDir = fullfile(extrDir, ['+', modelName], ['+', phaseName]);
            if ~exist(phaseDir, 'dir'), mkdir(phaseDir); end
            
            % Write inequality constraint function for this phase
            fid = fopen(fullfile(phaseDir, 'inequaliyConstraints.m'), 'w');
            fwrite(fid, strtrim(parts{i+1}));
            fclose(fid);
            
            % Copy class indices function from original model to maintain consistency
            origPhaseDir = fileparts(which([origModelName, '.', phaseName, '.get_classIndices']));
            if ~isempty(origPhaseDir)
                copyfile(fullfile(origPhaseDir, 'get_classIndices.m'), ...
                         fullfile(phaseDir, 'get_classIndices.m'));
            else
                warning('Could not find get_classIndices.m for phase %s in original model.', phaseName);
            end
        end
        addpath(extrDir);  % Add temporary directory to MATLAB path
    end
    
    % === Analyze each phase for active constraints ===
    for phaseIdx = 1:length(traj)
        phase = phaseSequence{phaseIdx};
        modelPhaseName = [modelName, '.', phase];
        
        % Convert trajectory to optimization variable format
        Z_i = interpolate_traj2Z(config, traj(phaseIdx), config.N(phaseIdx));
    
        % Load phase-specific functions
        ineq_func        = str2func([modelPhaseName, '.inequalityConstraints']);
        get_classIndices = str2func([modelPhaseName, '.get_classIndices']);
        
        % Get indices for states and controls in this phase
        [idx_x, idx_u] = get_classIndices();
    
        % Retrieve inequality constraint bounds and names
        [~, lb_all, ub_all, boundName_all] = ineq_func(config, Z_i, P_sym, config.N(phaseIdx));
    
        % Determine number of unique constraints (handles vectorized constraints)
        firstName = boundName_all{1};
        repeatStart = find(strcmp(boundName_all, firstName), 2, 'first');
        if length(repeatStart) == 2
            numUnique = repeatStart(2) - 1;  % Number of unique constraint types
        else
            numUnique = length(boundName_all);
        end
    
        % Extract unique constraint information
        boundName = boundName_all(1:numUnique);
        lb        = lb_all(1:numUnique);
        ub        = ub_all(1:numUnique);
    
        % Check each unique constraint for activity
        for cIdx = 1:numUnique
            name = boundName{cIdx};
    
            % --- Map constraint name to trajectory values ---
            stateMatch = find(strcmp(config.defaultStateNames, name), 1);
            inputMatch = find(strcmp(config.defaultInputNames, name), 1);
    
            if ~isempty(stateMatch)
                % State variable constraint
                idx = find(stateMatch == idx_x);
                vals = traj(phaseIdx).x(:, idx);
                % Store bound values for visualization
                trajLB(phaseIdx).x(:,idx) = lb(cIdx) * ones(length(vals), 1);
                trajUB(phaseIdx).x(:,idx) = ub(cIdx) * ones(length(vals), 1);
    
            elseif ~isempty(inputMatch)
                % Control input constraint
                idx = find(inputMatch == idx_u);
                vals = traj(phaseIdx).u(:, idx);
                trajLB(phaseIdx).u(:,idx) = lb(cIdx) * ones(length(vals), 1);
                trajUB(phaseIdx).u(:,idx) = ub(cIdx) * ones(length(vals), 1);
    
            elseif strcmp(name, 'dt')
                % Time step constraint
                N_phase     = config.N(phaseIdx);
                phaseDur    = traj(phaseIdx).t(end) - traj(phaseIdx).t(1);
                vals        = phaseDur * ones(1, length(traj(phaseIdx).t));
                lb(cIdx)    = lb(cIdx) * N_phase;  % Scale bounds for total phase duration
                ub(cIdx)    = ub(cIdx) * N_phase;
    
                trajLB(phaseIdx).t = lb(cIdx) * ones(length(vals),1);
                trajUB(phaseIdx).t = ub(cIdx) * ones(length(vals), 1);
    
            elseif any(strcmp(name, config.optParameterNames))
                % Optimized parameter constraint
                idx = find(strcmp(name, config.optParameterNames));
                vals = traj(phaseIdx).p(:,idx);
                trajLB(phaseIdx).p(:,idx) = lb(cIdx) * ones(length(vals), 1);
                trajUB(phaseIdx).p(:,idx) = ub(cIdx) * ones(length(vals), 1);
    
            elseif strcmp(name, 'ground_clearance')
                % Ground clearance constraint (y - l*cos(alpha))
                y_idx     = find(idx_x == find(strcmp(config.defaultStateNames, 'y')));
                alpha_idx = find(idx_x == find(strcmp(config.defaultStateNames, 'alpha')));
                l_idx     = find(idx_x == find(strcmp(config.defaultStateNames, 'l')));
                vals = traj(phaseIdx).x(:, y_idx) - ...
                       traj(phaseIdx).x(:, l_idx) .* cos(traj(phaseIdx).x(:, alpha_idx));
    
            elseif strcmp(name, 'jerk_alpha')
                % Hip torque jerk constraint
                u_alpha_idx = find(idx_u == find(strcmp(config.defaultInputNames, 'u_alpha')));
                dt = traj(phaseIdx).t(2) - traj(phaseIdx).t(1);
                vals = [0, diff(traj(phaseIdx).u(:, u_alpha_idx))']/dt;
    
            elseif strcmp(name, 'jerk_l')
                % Leg force jerk constraint
                u_l_idx = find(idx_u == find(strcmp(config.defaultInputNames, 'u_l')));
                dt = traj(phaseIdx).t(2) - traj(phaseIdx).t(1);
                vals = [0, diff(traj(phaseIdx).u(:, u_l_idx))']/dt;
                
            elseif strcmp(name, 'lambda')
                % Ground reaction force constraint
                vals = grf_func(traj(phaseIdx).x', traj(phaseIdx).u', P_num);
    
            else
                warning('Unknown bound name: %s', name);
                continue;
            end
    
            % --- Check if constraint is active (touching bounds) ---
            isActive = (vals <= lb(cIdx) + tol) | (vals >= ub(cIdx) - tol);
            if any(isActive)
                activeBounds(end+1) = struct( ...
                    'phase', phase, ...*
                    'name', name, ...
                    'idx', find(isActive) ...
                );
            end
        end
        
        % Initialize missing bound fields with infinite bounds
        if ~isfield(trajLB(phaseIdx), 'u') || isempty(trajLB(phaseIdx).u)
            trajLB(phaseIdx).u = -inf * ones(size(trajLB(phaseIdx).x, 1), numel(idx_u));
        end
        
        if ~isfield(trajUB(phaseIdx), 'u') || isempty(trajUB(phaseIdx).u)
            trajUB(phaseIdx).u = +inf * ones(size(trajUB(phaseIdx).x, 1), numel(idx_u));
        end
    end
    
    % === Clean up temporary files if external constraints were used ===
    if useExtr
        rmpath(extrDir);
        % Remove temporary directory and files
        if exist(extrDir, 'dir')
            rmdir(extrDir, 's');
        end
    end
end