function [eventFcn, jumpMap, nx_post] = createEventJump()
% CREATEEVENTJUMP Create CasADi functions for stance-to-stance transition
%
%   Constructs functions that define the nadir event (minimum leg length)
%   and the corresponding state jump map for stance-to-stance transitions.
%   Used for detecting the lowest point in the hopping cycle.
%
%   Inputs: None (uses internal model definitions)
%
%   Outputs:
%     eventFcn - CasADi function for nadir event detection
%     jumpMap  - CasADi function for state transition mapping
%     nx_post  - Dimension of post-transition state vector
%
%   Author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*
    import HF_noPitch_min.Stance.*;  % Import stance phase functions

    % Get model dimensions
    [nx, ~, nu, np] = get_sizes();

    % Define symbolic variables
    states = MX.sym('states', nx);      % Current state vector
    inputs = MX.sym('inputs', nu);      % Control inputs  
    parameters = MX.sym('parameters', np);  % Model parameters

    % Load pre-defined event and jump functions
    e = S2S.event(states, inputs, parameters);  % Nadir event: -dl (minimum leg length)
    g = S2S.jump(states, inputs, parameters);   % Jump map: continuous state [q; dq]

    nx_post = numel(g);  % Dimension of post-transition state

    % Create CasADi functions
    eventFcn = Function('event', {states, inputs, parameters}, {e});
    jumpMap = Function('jumpMap', {states, inputs, parameters}, {g});
end