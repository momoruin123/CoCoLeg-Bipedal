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
    import BF_min.SingleStance.*
    
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
    % Only one phase here (i = 1:1)
    for i = 1:nPhases
        phase = config.phaseSequence{i};
        N     = config.N(i);
        Z_i   = Z{i};

        % Get phase dimensions
        get_sizes = str2func([modelName, '.', phase, '.get_sizes']);
        [nx, ~, nu, ~] = get_sizes();

        % Extract the initial state and final state
        x_init = Z_i(1:nx);
        x_end  = Z_i(end-nx-nu-na+1:end-nu-na);

        % Extract time step if optimize time
        if config.optimizeTimeFlag
            dt = Z_i(nx+nu+1);
            % Accumulate total stride time
            totalTime = totalTime + N * dt;
        end
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
            % ----------------------------------------------------
            % NEED DEFINED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % ----------------------------------------------------
            % velocity constraints
            q0 = x_init(1:5);
            qT = x_end(1:5);

            v_target = config.operatingCond.v_avg;
            linkPts_init = LinkPositions_f([0;0;q0], params);
            linkPts_end = LinkPositions_f([0;0;qT], params);
            hipPt_init = linkPts_init(:,3);
            hipPt_end = linkPts_end(:,3);
            dx_hip = hipPt_end(1) - hipPt_init(1);
            eqs = dx_hip - v_target*totalTime;
            r = vertcat(r, eqs);

            % % Compute COM
            % % Using average hip joint average velocity as an indicator
            % LinkPos_init   = LinkPositions(x_init, params);
            % LinkPos_end    = LinkPositions(x_end, params);
            % hipPos_init_x  = LinkPos_init(1,3);
            % hipPos_end_x   = LinkPos_end(1,3);
            % eqs = hipPos_end_x - hipPos_init_x - config.operatingCond.v_avg * totalTime;
            % r = vertcat(r, eqs);

        otherwise
            error('No such optGaitMethod')
    end

    % Name the constraints
    n_r = size(r,1);
    for i = 1:n_r
        constraintName{i,1} = 'operating point';
    end        
end