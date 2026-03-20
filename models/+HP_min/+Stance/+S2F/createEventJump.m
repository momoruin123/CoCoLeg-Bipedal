function [eventFcn, jumpMap, nx_post] = createEventJump()
%CREATEEVENTJUMP Creates stance-to-flight transition functions
%    
%    author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*
    import HP_min.Stance.*;

    % Get system dimensions
    [nx, nq, nu, np] = get_sizes();

    % Symbolic variables
    states     = MX.sym('states', nx);
    inputs     = MX.sym('inputs', nu);
    parameters = MX.sym('parameters', np);

    % Get event and jump from S2F package
    e = S2F.event(states, inputs, parameters);
    g = S2F.jump(states, inputs, parameters);

    nx_post = numel(g);
    
    % Create CasADi functions
    eventFcn = Function('event', {states, inputs, parameters}, {e});
    jumpMap  = Function('jumpMap', {states, inputs, parameters}, {g});
end