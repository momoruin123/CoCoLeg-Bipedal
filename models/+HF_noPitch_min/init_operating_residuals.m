function [r, constraintName] = init_operating_residuals(config, Z)
% INIT_OPERATING_RESIDUALS Compute operating point and initial condition constraints
%
%   Enforces operating conditions for periodic locomotion including:
%   - Average velocity constraint over full stride
%   - Initial position constraint
%   - Foot rolling kinematics for stance phase transitions
%
%   Inputs:
%     config - system configuration structure with operating conditions
%     Z      - decision variables (cell array per phase or full vector)
%
%   Outputs:
%     r              - constraint residuals [operating_point; initial_conditions]
%     constraintName - names of each constraint element
%
%   author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025


    import casadi.*
    
    % Get model parameters and initialize variables
    params = getModelParameters(config, 0, 0);
    totalTime = 0;
    nPhases = numel(config.phaseSequence);
    modelName = config.model_name;
    na = numel(config.optParameterNames) + config.optimizeTimeFlag;
    
    % Convert monolithic Z vector to cell array per phase if needed
    if ~iscell(Z)
        Z_vec = Z;
        Z = {};
        idx_z = 1;
        for i = 1:nPhases
            phase = config.phaseSequence{i};
            N     = config.N(i);
            % Get phase dimensions
            get_sizes = str2func([modelName, '.', phase, '.get_sizes']);
            [nx, ~, nu, ~] = get_sizes();
            n_z = (N+1)*(nx+nu+na);
            Z{i} = Z_vec(idx_z: idx_z+n_z-1);
            idx_z = idx_z+n_z;
        end
    end
    
    % Extract key variables from each phase
    for i = 1:nPhases
        phase = config.phaseSequence{i};
        N     = config.N(i);
        Z_i   = Z{i};

        % Get phase dimensions
        get_sizes = str2func([modelName, '.', phase, '.get_sizes']);
        [nx, ~, nu, ~] = get_sizes();

        if i == 1
            % Extract initial state and time step from first phase
            x_init = Z_i(1:nx);
            dt     = Z_i(nx+nu+1);
        elseif i == nPhases
            % Extract final state and time step from last phase
            x_end  = Z_i(end-nx-nu-na+1:end-nu-na);
            dt     = Z_i(end-na+1);
        else 
            % Extract time step from intermediate phases
            dt = Z_i(nx+nu+1);
        end

        % Compute foot position at end of flight phase for rolling kinematics
        if strcmp(phase, 'Flight')
            x_1_end  = Z_i(end-nx-nu-na+1:end-nu-na);
            % Foot position: x + (phi_0 + phi)*r_f + l*sin(phi_0 + phi)
            d0 = x_1_end(1) + (params(getParameterIndex('phi_0', config)) + x_1_end(3)) * params(getParameterIndex('r_f', config)) + ...
                 x_1_end(4) * sin(params(getParameterIndex('phi_0', config)) + x_1_end(3));
        end

        % Accumulate total stride time
        totalTime = totalTime + N * dt;
    end
    
    % Initialize constraint arrays
    r = MX.zeros(5*nx,1);
    constraintName = cell(5*nx,1);
    idx_r = 0;
    
    % Operating point constraint: enforce average velocity
    % Compute final foot position using rolling contact kinematics
    x_star = d0 - params(getParameterIndex('r_f', config)) * (params(getParameterIndex('phi_0', config)) + x_end(1)) - ...
             x_end(2) * sin(params(getParameterIndex('phi_0', config)) + x_end(1));
    
    % Average velocity constraint: x_final - x_initial = v_avg * totalTime
    r(idx_r+1) = x_star - config.operatingCond.x0 - config.operatingCond.v_avg * totalTime;
    constraintName{idx_r+1} = 'operating point';
    idx_r = idx_r+1;
    
    % Initial position constraint
    r(idx_r+1) = x_init(1) - config.operatingCond.x0;  
    constraintName{idx_r+1} = 'initial conditions';
    idx_r = idx_r+1;
    
    % Truncate to actual number of constraints
    r = r(1:idx_r);
    constraintName = constraintName(1:idx_r);
end