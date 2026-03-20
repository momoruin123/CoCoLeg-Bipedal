function Dyn = createDynamics()
    import casadi.*
    import BF_run.Flight.*;

    % Get system dimensions
    [nx, nq, nu, np] = get_sizes();

    % Symbolic variables
    states     = MX.sym('states', nx);
    inputs     = MX.sym('inputs', nu);
    parameters = MX.sym('parameters', np);

    % Compute dynamics: M*ddq = forces
    h  = forces_f(states, inputs, parameters);
    M  = massMatrix_f(states, parameters);

    h = MX(h);
    M = MX(M);

    acc = M \ h;
    dq = states(nq+1:nx);
    dstate = [dq; acc];
    
    % Create CasADi function
    Dyn = Function('Dyn', {states, inputs, parameters}, {dstate});
end