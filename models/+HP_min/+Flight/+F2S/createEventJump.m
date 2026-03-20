function [eventFcn, jumpMap, nx_post] = createEventJump()
%CREATEEVENTJUMP Creates flight-to-stance transition functions
%    
%    author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*
    import HP_min.Flight.*;

    % Get system dimensions
    [nx, nq, nu, np] = get_sizes();

    % Symbolic variables
    states     = MX.sym('states', nx);
    inputs     = MX.sym('inputs', nu);
    parameters = MX.sym('parameters', np);

    % Get event and jump from F2S package
    e = F2S.event(states, inputs, parameters);
    g = F2S.jump(states, inputs, parameters);

    nx_post = numel(g);
    
    % Create CasADi functions
    eventFcn = Function('eventFcn', {states, inputs, parameters}, {e});
    jumpMap = Function('jumpMap', {states, inputs, parameters}, {g});
end