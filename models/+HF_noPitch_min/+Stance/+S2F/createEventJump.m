function [eventFcn, jumpMap, nx_post] = createEventJump()
% CREATEEVENTJUMP Create CasADi functions for phase transition events and jumps
%
%   Constructs functions that define the lift-off event (transition from 
%   stance to flight) and the corresponding state jump map for the 
%   minimal coordinate representation.
%
%   Inputs: None (uses internal model definitions)
%
%   Outputs:
%     eventFcn - CasADi function for lift-off event detection
%     jumpMap  - CasADi function for state transition mapping
%     nx_post  - Dimension of post-transition state vector
%
%   Author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*
    import HF_noPitch_min.Stance.*;  % Import stance phase functions

    % Get model dimensions
    [nx, nq, nu, np] = get_sizes();

    % Define symbolic variables
    states = MX.sym('states', nx);      % Current state vector
    inputs = MX.sym('inputs', nu);      % Control inputs
    parameters = MX.sym('parameters', np);  % Model parameters

    % Compute state derivatives (dynamics)
    h = forces(states, inputs, parameters);  % Generalized forces
    M = massMatrix(states, parameters);      % Mass matrix
    dq = states(nq+1:end);                   % Current velocities
    acc = M \ h;                             % Solve for accelerations
    f_S = [dq; acc(1:nq)];                   % Full state derivative

    % Event function: net vertical force becoming zero (lift-off condition)
    threeBodiesAccjac_Y = S2F.threeBodiesAccjac_Y(states, parameters);  % Acceleration Jacobian
    masses = S2F.masses(parameters);          % Body masses [m_t, m_l, m_f]
    e = masses * (threeBodiesAccjac_Y * f_S + parameters(2) * ones(3,1));  % Net vertical force

    % Jump map: state transformation to minimal coordinates (x-position ignored)
    g = S2F.jump_noX(states, inputs, parameters);  % Minimal state representation
    nx_post = numel(g);  % Dimension of post-transition state

    % Create CasADi functions
    eventFcn = Function('event', {states, inputs, parameters}, {e});
    jumpMap = Function('jumpMap', {states, inputs, parameters}, {g});
end