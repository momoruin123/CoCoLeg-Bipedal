function [r, constraintName] = transition_residuals(config, phase, phaseTo, Z_pre, Z_post, P, periodicityFlag)
% TRANSITION_RESIDUALS Compute transition constraints between phases
%
%   Handles both regular phase transitions and periodic boundary conditions
%   Enforces event conditions, jump maps, and parameter continuity
%
%   Inputs:
%     config          - system configuration structure
%     phase           - current phase name
%     phaseTo         - next phase name  
%     Z_pre           - decision variables at end of current phase
%     Z_post          - decision variables at start of next phase
%     P               - model parameters
%     periodicityFlag - flag for periodic boundary conditions
%
%   Outputs:
%     r               - transition constraint residuals
%     constraintName  - names of constraint elements
%
%   author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*
    
    % Calculate number of auxiliary variables (time + optimized parameters)
    na = numel(config.optParameterNames) + config.optimizeTimeFlag;
    
    % Construct package names for current phase and transition
    packageName = [config.model_name, '.', phase];
    transPackageName = [config.model_name, '.', phase, '.', phase(1), '2', phaseTo(1)];
    
    % Initialize outputs
    r = MX.zeros(numel(Z_pre) + numel(Z_post));
    constraintName = cell(numel(Z_pre) + numel(Z_post), 1);
    idx_r = 0;  % Residual index counter
    
    % Load model-specific functions
    get_sizes = str2func([packageName, '.get_sizes']);
    createEventJump = str2func([transPackageName, '.createEventJump']);
    
    % Get model dimensions and transition functions
    [nx, ~, nu, ~] = get_sizes();
    [eventFcn, jumpMap, nx_plus] = createEventJump();
    
    % Set selection matrix for periodicity (identity for regular transitions)
    if periodicityFlag
        S = config.periodicityMatrix;  % Periodicity matrix for cyclic motions
    else
        S = eye(nx_plus);       % Identity matrix for regular transitions
    end
    
    % Extract variables at the end of current phase (pre-transition)
    Z_pre = Z_pre(end - (nx + nu + na - 1):end);  % Last column of decision variables
    xk_eminus = Z_pre(1:nx, end);                 % State before transition
    uk_eminus = Z_pre(nx + 1:nx + nu, end);       % Control before transition  
    ak_eminus = Z_pre(nx + nu + 1:nx + nu + na, end);  % Auxiliary variables before transition

    % Construct full parameter vector (fixed + optimized parameters)
    p_aux = ak_eminus(end - na + 2:end);  % Optimized parameters
    p_full = [P; p_aux];                  % Combined parameter vector

    % Event function constraint - ensures transition condition is met
    r(idx_r + 1) = eventFcn(xk_eminus, uk_eminus, p_full);
    constraintName{idx_r + 1} = 'event';
    idx_r = idx_r + 1;

    % Extract variables at the start of next phase (post-transition)
    xk_eplus = Z_post(1:nx_plus);                     % State after transition
    ak_eplus = Z_post(nx_plus + nu + 1:nx_plus + nu + na);  % Auxiliary variables after transition

    % Jump map constraint - ensures state continuity through transition
    r(idx_r + (1:size(S, 1))) = S * (jumpMap(xk_eminus, uk_eminus, p_full) - xk_eplus);
    constraintName(idx_r + (1:size(S, 1))) = {'jump map'};
    idx_r = idx_r + size(S, 1);

    % Parameter continuity constraint (only for non-periodic transitions)
    if ~periodicityFlag
        % Ensure optimized parameters remain constant across phase boundaries
        r_param = ak_eplus(end - na + 2:end) - ak_eminus(end - na + 2:end);
        idx = idx_r + 1:idx_r + numel(r_param);
        r(idx(:)) = reshape(r_param, size(r(idx(:))));
        constraintName(idx_r + 1:idx_r + numel(r_param)) = {'constant parameter'};
        idx_r = idx_r + numel(r_param);
    end

    % Truncate outputs to actual sizes
    r = r(1:idx_r);
    constraintName = constraintName(1:idx_r);
end