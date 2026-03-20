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
    import BF_run.SingleStance.*
    
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

        % Accumulate total stride time
        totalTime = totalTime + N * dt;
    end

    % Initialize constraint arrays
    r = MX.zeros(0,1);
    constraintName = cell(0,1);

    % Emplement Gait optimizing methods
    switch config.optGaitMethods
        % optimize gait with fix step length
        case 'Fix_step_length'
            L = config.stepLength;
            % End position of swing foot (2x1)
            SwingFoot = cons_SwingFoot(x_end, params);
            eqs = [SwingFoot(1) - L; SwingFoot(2) - 0];
            r = vertcat(r, eqs);

        case 'Fix_average_velocity'
            % velocity constraints
            q0 = vertcat(0,0,x_init(1:5));
            qT = x_end(1:7);

            v_target = config.operatingCond.v_avg;
            linkPts_init = LinkPositions_f(q0, params);
            linkPts_end = LinkPositions_f(qT, params);
            hipPt_init = linkPts_init(:,3);
            hipPt_end = linkPts_end(:,3);
            dx_hip = hipPt_end(1) - hipPt_init(1);
            eqs = dx_hip - v_target*totalTime;
            r = vertcat(r, eqs);

        otherwise
            error('No such optGaitMethod')
    end

    % Name the constraints
    n_r = size(r,1);
    for i = 1:n_r
        constraintName{i,1} = 'operating point';
    end        
end