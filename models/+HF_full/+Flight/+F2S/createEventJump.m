function [eventFcn, jumpMap, nx_post] = createEventJump()
%CREATEEVENTJUMP Creates flight-to-stance transition functions

    import casadi.*
    import HF_full.Flight.*;

    % Get system dimensions
    [nx, nq, nu, np] = get_sizes();

    % Symbolic variables
    states     = MX.sym('states', nx);
    inputs     = MX.sym('inputs', nu);
    parameters = MX.sym('parameters', np);

    % Compute impact dynamics
    M  = massMatrix(states, parameters);            % Mass matrix
    W  = contactJacobian(states, parameters);       % Contact Jacobian

    % Solve impact equation: [M -W; W' 0] * [dq_plus; impulse] = [M*dq; 0]
    q  = states(1:nq);
    dq = states(nq+1:end);
    plus = [M -W; W' zeros(size(W,2))] \ [M*dq; zeros(size(W,2),1)];

    % Event and jump map
    e = F2S.event(states, inputs, parameters);  % Touchdown event
    g = [q; plus(1:nq)];                       % Post-impact state (positions unchanged)

    nx_post = numel(g);
    
    % Create CasADi functions
    eventFcn = Function('eventFcn', {states, inputs, parameters}, {e});
    jumpMap = Function('jumpMap', {states, inputs, parameters}, {g});
end