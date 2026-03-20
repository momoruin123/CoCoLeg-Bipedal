function [Z_sym, P, w_all, r_all, constraintNames, h_all, bounds_all, boundNames] = getConstraintsAndCosts(config, traj)
% GETCONSTRAINTSANDCOSTS Generate symbolic constraints and cost functions for optimization
%
%   Inputs:
%     config - system configuration structure
%     traj   - trajectory data for initialization
%
%   Outputs:
%     Z_sym          - symbolic decision variables
%     P              - symbolic model parameters
%     w_all          - total cost function
%     r_all          - equality constraints (dynamics + transitions)
%     constraintNames- names of equality constraints
%     h_all          - inequality constraints
%     bounds_all     - inequality constraint bounds [lb, ub]
%     boundNames     - names of inequality constraints
%  
%   author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*

    % Extract configuration parameters
    modelName     = config.model_name;
    phaseSequence = config.phaseSequence;
    useIneq       = config.useInequalityConstraints;
    N_array       = config.N;
    
    % Initialize function handles
    init_operating_residuals = str2func([modelName, '.init_operating_residuals']);

    % Get symbolic parameters
    P = getModelParameters(config, 1, 1);

    % Initialize symbolic decision variables for each phase
    nPhases = numel(phaseSequence);
    Z_sym = cell(1, nPhases);
    Z0 = cell(1, nPhases);
    nz = 0;

    for i = 1:nPhases
        Z0{i} = interpolate_traj2Z(config, traj(i), N_array(i));
        Z_sym{i} = MX.sym(['Z_', phaseSequence{i}, '_', num2str(i)], numel(Z0{i}));
        nz = nz + numel(Z_sym{i});
    end

    % Initialize constraint and cost arrays
    r_all = MX.zeros(nz);
    h_all = MX.zeros(2*nz);
    w_all = MX.zeros(1);
    constraintNames = cell(nz, 1);
    lb_all = cell(1, nPhases);
    ub_all = cell(1, nPhases);
    boundNames = cell(2*nz, 1);

    % Initialize constraint counters
    idx_r = 0;
    idx_h = 0;

    % Add operating point residuals
    [r_init, name_init] = init_operating_residuals(config, Z_sym);
    r_all(idx_r+1:idx_r+numel(r_init)) = r_init;
    constraintNames(idx_r+1:idx_r+numel(r_init)) = name_init;
    idx_r = idx_r + numel(r_init);
    
    % Process each phase
    for i = 1:nPhases
        % Update inequality constraint coefficients if specified
        if isfield(config, 'ineqConstraints') && isfield(config.ineqConstraints, 'coef') && isfield(config.ineqConstraints, 'array') 
           config.ineqConstraints.coef = config.ineqConstraints.array{i};
        end
        
        phase = phaseSequence{i};
        Z_i = Z_sym{i};
        modelPhaseName = [modelName, '.', phase];

        % Add dynamics residuals and cost
        [r_dyn, w_dyn, name_dyn] = dynamics_residuals(config, modelPhaseName, Z_i, P, N_array(i));
        r_all(idx_r+1:idx_r+numel(r_dyn)) = r_dyn;
        name_dyn = cellfun(@(s) [s ' ' phase], name_dyn, 'UniformOutput', false);
        constraintNames(idx_r+1:idx_r+numel(r_dyn)) = name_dyn;
        w_all = w_all + w_dyn;
        idx_r = idx_r + numel(r_dyn);

        % Add inequality constraints if enabled
        if useIneq
            try
                ineq_func = str2func([modelPhaseName, '.inequalityConstraints']);
                [h_i, lb_i, ub_i, name_i] = ineq_func(config, Z_i, P, N_array(i));
                
                h_all(idx_h+1:idx_h+numel(h_i)) = h_i;
                boundNames(idx_h+1:idx_h+numel(h_i)) = name_i;
                idx_h = idx_h + numel(h_i);
            catch ME
                warning('Failed to evaluate inequality constraints for phase "%s": %s', ...
                    modelPhaseName, ME.message);
                rethrow(ME);
            end
        else
            h_i = [];
            lb_i = [];
            ub_i = [];
            name_i = {};
        end

        lb_all{i} = lb_i;
        ub_all{i} = ub_i;

        % Add transition residuals
        if i < nPhases
            % Transition to next phase
            phaseTo = phaseSequence{i+1};
            Z_next = Z_sym{i+1};
            periodicityFlag = 0;
        else
            % Periodic transition back to first phase
            phaseTo = phaseSequence{1};
            Z_next = Z_sym{1};
            periodicityFlag = 1 & config.periodicityFlag;
        end
        [r_trans, name_trans] = transition_residuals(config, phase, phaseTo, Z_i, Z_next, P, periodicityFlag);  
        r_all(idx_r+1:idx_r+numel(r_trans)) = r_trans;
        constraintNames(idx_r+1:idx_r+numel(r_trans)) = name_trans;
        idx_r = idx_r + numel(r_trans);
    end

    % Concatenate bounds and normalize cost
    lb_all = vertcat(lb_all{:});
    ub_all = vertcat(ub_all{:});
    bounds_all = [lb_all, ub_all];

    w_all = normalize_cost(config, w_all, Z_sym, P);
    
    % Truncate arrays to actual sizes
    r_all = r_all(1:idx_r);
    constraintNames = constraintNames(1:idx_r);
    h_all = h_all(1:idx_h);
    boundNames = boundNames(1:idx_h);
    Z_sym = vertcat(Z_sym{:});
end
