function [Z_final, cost, return_status] = continuation_single_step(Z_init, grid_parfor, traj, config)
% CONTINUATION_SINGLE_STEP Perform single optimization step in parameter continuation
%
%   Executes one optimization step for a specific grid point in a parameter
%   continuation study. Updates operating conditions and parameters, then
%   solves the optimization problem using warm-start initialization.
%
%   Inputs:
%     Z_init      - initial guess for decision variables
%     grid_parfor - grid point values [operating conditions, parameters]
%     traj        - trajectory structure for constraint formulation
%     config      - system configuration structure
%
%   Outputs:
%     Z_final     - optimized decision variables
%     cost        - final cost value
%     return_status - IPOPT solver return status
%
%   author: C. David Remy, Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    % Safe defaults so caller can continue even if this step fails.
    Z_final = full(Z_init(:)).';
    cost = inf;
    return_status = "Failed_Runtime";

    % Set model stiffness parameters
    p_vec = getModelParameters(config, 0, 0);
    sf = p_vec(1)*p_vec(2)*p_vec(3);

    try
        % Update operating conditions from grid point values
        n_opCond = numel(config.cont.gridOpCondNames);
        n_param = numel(config.cont.gridParamNames);
        
        % Set operating condition values
        for i = 1:n_opCond
            opCondName = config.cont.gridOpCondNames{i};
            config.operatingCond.(opCondName) = grid_parfor(i);
        end
        
        % Set parameter values  
        for i = 1:n_param
            paramName = config.cont.gridParamNames{i};
            config.paramValues.(paramName) = grid_parfor(n_opCond + i)/sf;
        end
        
        % Generate optimization problem constraints and costs
        [Z, P, w, r, ~, h, b, ~] = getConstraintsAndCosts(config, traj);
        
        % Solve optimization problem
        [solutionIPOPT, stats] = solveCasadi(config, Z, P, r, h, w, b, Z_init);
        
        % Extract results
        return_status = string(stats.return_status);
        Z_final = full(solutionIPOPT.x).';
        cost = full(solutionIPOPT.f);
    catch
        % Keep fallback outputs to mark this grid point as failed.
    end
end
