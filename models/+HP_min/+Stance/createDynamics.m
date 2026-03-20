function Dyn = createDynamics()
%CREATEDYNAMICS Creates CasADi function for system dynamics
% Symbolic implementation for efficient optimization.
% Uses CasADi for automatic differentiation - provides exact derivatives 
% Avoids differentiating complex numerical MATLAB functions.
% Essential for reliable convergence in trajectory optimization.
%
%   Author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025


    import casadi.*
    import HP_min.Stance.*;

    % Get system dimensions
    [nx, nq, nu, np] = get_sizes();

    % Symbolic variables
    states     = MX.sym('states', nx);
    inputs     = MX.sym('inputs', nu);
    parameters = MX.sym('parameters', np);

    % Compute dynamics: M*ddq = forces
    h  = forces(states, inputs, parameters);
    M  = massMatrix(states, parameters);
    acc = M\h;
    dq = states(nq+1:nx);
    dstate = [dq; acc];
    
    % Create CasADi function
    Dyn = Function('Dyn', {states, inputs, parameters}, {dstate});
end