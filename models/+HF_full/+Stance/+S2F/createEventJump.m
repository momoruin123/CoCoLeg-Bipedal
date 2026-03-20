function [eventFcn, jumpMap, nx_post] = createEventJump()
%CREATEEVENTJUMP Creates event function and jump map for phase transitions

    import casadi.*
    import HF_full.Stance.*;

    % Get system dimensions
    [nx, nq, nu, np] = get_sizes();

    % Symbolic variables
    states     = MX.sym('states', nx);
    inputs     = MX.sym('inputs', nu);
    parameters = MX.sym('parameters', np);

    % Compute constrained dynamics
    h  = forces(states, inputs, parameters);        % Generalized forces
    M  = massMatrix(states, parameters);            % Mass matrix
    W  = contactJacobian(states, parameters);       % Contact Jacobian
    dW = contactJacobianDerivative(states, parameters); % Contact Jacobian derivative
    
    % Solve constrained dynamics: [M -W; W' 0] * [ddq; lambda] = [h; -dW'*dq]
    dq  = states(nq+1:end);
    acc = [M -W; W' zeros(size(W,2))]\[h; -dW'*dq];
    
    % Event: normal contact force (last Lagrange multiplier)
    event = acc(end);

    % Jump map from S2F package
    g = S2F.jump(states, inputs, parameters);
    nx_post = numel(g);
    
    % Create CasADi functions
    eventFcn = Function('event', {states, inputs, parameters}, {event});
    jumpMap  = Function('jumpMap', {states, inputs, parameters}, {g});
end