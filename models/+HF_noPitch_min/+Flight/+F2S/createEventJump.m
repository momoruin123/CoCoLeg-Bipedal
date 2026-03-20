function [eventFcn, jumpMap, nx_post] = createEventJump()
% CREATEEVENTJUMP Create CasADi functions for flight-to-stance transition
%
%   Constructs functions that define the touchdown event and the corresponding
%   state jump map for flight-to-stance transitions. Handles the impact
%   dynamics and coordinate reduction when the foot contacts the ground.
%
%   Inputs: None (uses internal model definitions)
%
%   Outputs:
%     eventFcn - CasADi function for touchdown event detection
%     jumpMap  - CasADi function for state transition mapping
%     nx_post  - Dimension of post-transition state vector
%
%   Author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*
    import HF_noPitch_min.Flight.*;  % Import flight phase functions

    % Get model dimensions
    [nx, nq, nu, np] = get_sizes();

    % Define symbolic variables
    states = MX.sym('states', nx);      % Current state vector
    inputs = MX.sym('inputs', nu);      % Control inputs
    parameters = MX.sym('parameters', np);  % Model parameters

    % Compute dynamics components
    h = forces(states, inputs, parameters);  % Generalized forces
    M = massMatrix(states, parameters);      % Mass matrix
    W = contactJacobian(states, parameters); % Contact Jacobian
    dW = contactJacobianDerivative(states, parameters); % Contact Jacobian derivative

    % Extract positions and velocities
    q = states(1:nq);    % Generalized coordinates
    dq = states(nq+1:end); % Generalized velocities

    % Solve impact dynamics: [M, -W; W', 0] * [dq_plus; lambda] = [M*dq_minus; 0]
    % for post-impact velocities and contact forces
    plus = [M, -W; W', zeros(size(W,2))] \ [M*dq; zeros(size(W,2),1)];

    % Event and jump functions
    e = F2S.event(states, inputs, parameters);  % Touchdown event (foot height = 0)
    g = [q(3:nq); plus(3:nq)];  % Post-impact state: [alpha, l, dalpha, dl] (minimal coordinates)

    nx_post = numel(g);  % Dimension of post-transition state

    % Create CasADi functions
    eventFcn = Function('eventFcn', {states, inputs, parameters}, {e});
    jumpMap = Function('jumpMap', {states, inputs, parameters}, {g});
end