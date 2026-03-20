function [eventFcn, jumpMap, nx_post] = createEventJump()
    import casadi.*
    % Import all public members from mypackage
    import HF_full.Flight.*;

    % find a way to get nx, nq, nu and np.
    [nx, nq, nu, np] = get_sizes();

    states     = MX.sym('states', nx);
    inputs     = MX.sym('inputs', nu);
    parameters = MX.sym('parameters', np);

    e = F2F.event(states,inputs,parameters);
    g = F2F.jump(states, inputs, parameters);

    nx_post = numel(g);
    
    eventFcn = Function('eventFcn', {states, inputs, parameters}, {e});
    jumpMap = Function('jumpMap', {states, inputs, parameters}, {g});
end