function Dyn = createDynamics()
%CREATEDYNAMICS Creates CasADi function for system dynamics
% Symbolic implementation for efficient optimization.
% Uses CasADi for automatic differentiation - provides exact derivatives 
% Avoids differentiating complex numerical MATLAB functions.
% Essential for reliable convergence in trajectory optimization.
%
%   Author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*
    import HF_noPitch_min.Stance.*;  % Import stance phase functions

    % Get model dimensions
    [nx, nq, nu, np] = get_sizes();

    % Define symbolic variables
    states = MX.sym('states', nx);      % State vector [q; dq]
    inputs = MX.sym('inputs', nu);      % Control input vector
    parameters = MX.sym('parameters', np);  % Parameter vector

    % Compute forces and mass matrix
    h = forces(states, inputs, parameters);  % Generalized forces
    M = massMatrix(states, parameters);      % Mass matrix

    % Extract velocities and compute accelerations
    dq = states(nq+1:end);           % Generalized velocities
    acc = M \ h;                     % Solve M*ddq = h for accelerations
    dstate = [dq; acc(1:nq)];        % State derivative [dq; ddq]

    % Create CasADi function
    Dyn = Function('Dyn', {states, inputs, parameters}, {dstate});
end