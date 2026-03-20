function [eventFcn, jumpMap, nx_post] = createEventJump()
% CREATEEVENTJUMP Create CasADi functions for flight-to-flight transition
%
%   Constructs functions that define the apex event and the corresponding
%   state jump map for flight-to-flight transitions. Used for detecting
%   the highest point in the flight phase for phase splitting.
%
%   Inputs: None (uses internal model definitions)
%
%   Outputs:
%     eventFcn - CasADi function for apex event detection
%     jumpMap  - CasADi function for state transition mapping
%     nx_post  - Dimension of post-transition state vector

    import casadi.*
    import HF_noPitch_min.Flight.*;  % Import flight phase functions

    % Get model dimensions
    [nx, nq, nu, np] = get_sizes();

    % Define symbolic variables
    states = MX.sym('states', nx);      % Current state vector
    inputs = MX.sym('inputs', nu);      % Control inputs
    parameters = MX.sym('parameters', np);  % Model parameters

    % Load pre-defined event and jump functions
    e = F2F.event(states, inputs, parameters);  % Apex event: -dy (maximum height)
    g = F2F.jump(states, inputs, parameters);   % Jump map: continuous state

    nx_post = numel(g);  % Dimension of post-transition state

    % Create CasADi functions
    eventFcn = Function('eventFcn', {states, inputs, parameters}, {e});
    jumpMap = Function('jumpMap', {states, inputs, parameters}, {g});
end