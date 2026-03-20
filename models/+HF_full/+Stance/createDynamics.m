function Dyn = createDynamics()
%CREATEDYNAMICS Creates CasADi function for system dynamics
% Symbolic implementation for efficient optimization.
% Uses CasADi for automatic differentiation - provides exact derivatives 
% Avoids differentiating complex numerical MATLAB functions.
% Essential for reliable convergence in trajectory optimization.
%
%   Author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*
    import HF_full.Stance.*;

    % Get system dimensions
    [nx, nq, nu, np] = get_sizes();

    % Symbolic variables
    states     = MX.sym('states', nx);      % [q; dq]
    inputs     = MX.sym('inputs', nu);      % controls
    parameters = MX.sym('parameters', np);  % physical parameters
    dq = states(nq+1:nx);                   % Extract velocities

    % Compute dynamics: M*ddq = forces
    h  = forces(states, inputs, parameters);  % Generalized forces
    M  = massMatrix(states, parameters);      % Mass matrix
    W  = contactJacobian(states, parameters);       % Contact Jacobian
    dW = contactJacobianDerivative(states, parameters); % Contact Jacobian derivative
    acc_lamb = [M -W; W' zeros(size(W,2))] \ [h;-dW'*dq];
    
    % State derivative: [velocities; accelerations]
    dstate = [dq; acc_lamb(1:nq)];                     % Full state derivative
    
    % Create CasADi function
    Dyn = Function('Dyn', {states, inputs, parameters}, {dstate});
end