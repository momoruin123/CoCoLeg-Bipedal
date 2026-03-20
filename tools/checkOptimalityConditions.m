function [jacobians_Sym, optimalityConditions_Num] = checkOptimalityConditions(config, solutionIPOPT, Z, P, w, r, h)
% CHECKOPTIMALITYCONDITIONS Verify Karush-Kuhn-Tucker (KKT) optimality conditions
%
%   Analyzes the solution of a nonlinear optimization problem by computing
%   and evaluating the Karush-Kuhn-Tucker (KKT) optimality conditions.
%   Provides symbolic Jacobians and numerical validation of first-order
%   and second-order optimality conditions.
%
%   Inputs:
%     config         - system configuration structure
%     solutionIPOPT  - IPOPT solution structure containing optimal variables
%     Z              - symbolic decision variables
%     P              - symbolic model parameters
%     w              - cost function expression
%     r              - equality constraint expressions (dynamics + transitions)
%     h              - inequality constraint expressions
%
%   Outputs:
%     jacobians_Sym          - structure containing symbolic Jacobian functions:
%                              dr_dZ_eval: Jacobian of equality constraints
%                              dh_dZ_eval: Jacobian of inequality constraints
%                              de_dZ_eval: Gradient of cost function
%                              dL_dZ_eval: Gradient of Lagrangian
%                              ddL_dZ_eval: Hessian of Lagrangian
%     optimalityConditions_Num - structure containing numerical optimality analysis:
%                              dr_dZ, dh_dZ: Constraint Jacobians evaluated at solution
%                              dL_dZ: Lagrangian gradient (should be near zero)
%                              equality_lambda: Lagrange multipliers for equality constraints
%                              inequality_lambda: Lagrange multipliers for inequality constraints
%                              violatedEqualityConstIdx: Indices of violated equality constraints
%                              activeInequalityConstIdx: Indices of active inequality constraints
%                              ddL_dZ: Hessian of Lagrangian
%                              projectedJacobian: Reduced Hessian in null space of constraints
%                              minEigSecondOrder: Minimum eigenvalue of reduced Hessian
%                              EigSecondOrder: All eigenvalues of reduced Hessian
%
%  author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*

    % Get numerical parameter values
    p_num = getModelParameters(config, 0, 1);

    % Combine all constraints
    g = [r; h];
    
    % Define Lagrange multipliers for the Lagrangian
    lag_lambda = MX.sym('lag_lambda', numel(g));
    L = w + lag_lambda' * g;
    
    % Initialize symbolic Jacobians structure
    jacobians_Sym = {};
    
    %% Compute constraint Jacobians
    % Equality constraint Jacobian: dr/dZ
    dr_dZ = jacobian(r, Z);
    jacobians_Sym.dr_dZ_eval = casadi.Function('dr_dZ_eval', {Z, P}, {dr_dZ});
    
    % Inequality constraint Jacobian: dh/dZ
    dh_dZ = jacobian(h, Z);
    jacobians_Sym.dh_dZ_eval = casadi.Function('dh_dZ_eval', {Z, P}, {dh_dZ});
    
    %% Compute cost function gradient: dw/dZ
    de_dZ = jacobian(w, Z);
    jacobians_Sym.de_dZ_eval = casadi.Function('de_dZ_eval', {Z, P}, {de_dZ});
 
    %% Compute gradient of the Lagrangian: dL/dZ
    dL_dZ = jacobian(L, Z);
    jacobians_Sym.dL_dZ_eval = casadi.Function('dL_dZ_eval', {Z, P, lag_lambda}, {dL_dZ});

    %% Compute Hessian of the Lagrangian: d²L/dZ²
    ddL_dZ = hessian(L, Z);
    jacobians_Sym.ddL_dZ_eval = casadi.Function('ddL_dZ_eval', {Z, P, lag_lambda}, {ddL_dZ});

    %% Evaluate optimality conditions numerically at the solution
    optimalityConditions_Num = {};
    
    % Evaluate constraint Jacobians at optimal solution
    optimalityConditions_Num.dr_dZ = full(jacobians_Sym.dr_dZ_eval(solutionIPOPT.x, p_num));
    optimalityConditions_Num.dh_dZ = full(jacobians_Sym.dh_dZ_eval(solutionIPOPT.x, p_num));

    % Evaluate Lagrangian gradient at optimal solution
    optimalityConditions_Num.dL_dZ = full(jacobians_Sym.dL_dZ_eval(solutionIPOPT.x, p_num, solutionIPOPT.lam_g));

    % Analyze constraint satisfaction and activity
    optimalityConditions_Num.equality_lambda = full(solutionIPOPT.g);
    optimalityConditions_Num.violatedEqualityConstIdx = find(abs(optimalityConditions_Num.equality_lambda(1:size(optimalityConditions_Num.dr_dZ,1))) > 1e-9);
    
    optimalityConditions_Num.inequality_lambda = full(solutionIPOPT.lam_g);
    optimalityConditions_Num.activeInequalityConstIdx = find(abs(optimalityConditions_Num.inequality_lambda(numel(r)+1:end)) > 1e-4);

    % Second-order optimality conditions
    optimalityConditions_Num.ddL_dZ = full(jacobians_Sym.ddL_dZ_eval(solutionIPOPT.x, p_num, solutionIPOPT.lam_g));

    % Project Hessian onto null space of active constraints for reduced space analysis
    null_jacob = computeNull(optimalityConditions_Num.dr_dZ);
    res = null_jacob' * optimalityConditions_Num.ddL_dZ * null_jacob;
    
    optimalityConditions_Num.projectedJacobian = res;
    optimalityConditions_Num.minEigSecondOrder = min(eig(res));
    optimalityConditions_Num.EigSecondOrder = eig(res);
end

function t = computeNull(A)
% COMPUTENULL Compute null space basis using QR decomposition
%
%   Computes an orthonormal basis for the null space of matrix A using
%   economy-size QR decomposition. Assumes A has full row rank.
%
%   Input:
%     A - matrix to compute null space for (assumed full row rank)
%
%   Output:
%     t - orthonormal basis for the null space of A
%
% author: Maximillian Raff

    % Determine dimension of null space
    dim = size(A,2) - size(A,1);
    
    try
        % Compute QR decomposition of transpose for null space calculation
        [Q,~] = qr(sparse(A'));
    catch 
        % Enter debug mode if QR decomposition fails
        keyboard; % Pauses execution for interactive debugging
    end
    
    % Extract null space basis from last 'dim' columns of Q
    t = Q(:,(end+1-dim):end);
end
