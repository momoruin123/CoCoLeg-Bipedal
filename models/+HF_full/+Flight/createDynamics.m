function Dyn = createDynamics()

%CREATEDYNAMICS Creates CasADi function for system dynamics
% Symbolic implementation for efficient optimization.
% Uses CasADi for automatic differentiation - provides exact derivatives 
% Avoids differentiating complex numerical MATLAB functions.
% Essential for reliable convergence in trajectory optimization.
%
%   Author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*
    % Import all public members from mypackage
    import HF_full.Flight.*;

    [nx, nq, nu, np] = get_sizes();

    states     = MX.sym('states', nx);
    inputs     = MX.sym('inputs', nu);
    parameters = MX.sym('parameters', np);

    % [parameter, ~] = updateParameters(decisionParam, optParam);
    h  = forces(states, inputs, parameters);
    M  = massMatrix(states, parameters);                    % M = massMatrix(x, p_full);
    
    acc = M\h;
    dq = states(nq+1:nx);

    dstate = [dq; acc]; 
    
    Dyn = Function('Dyn', {states, inputs, parameters}, {dstate});
end