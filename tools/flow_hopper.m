function traj = flow_hopper(config, dt)
% FLOW Simulate hopper trajectory through multiple phases
%   traj = flow_hopper(config, dt) simulates the system dynamics through the 
%   specified phase sequence using the given configuration and time step
%
%   Inputs:
%     config - system configuration structure
%     dt     - time step or vector of time steps for each phase
%
%   Output:
%     traj   - trajectory structure containing states, controls, time, 
%              parameters, and constraint forces
%    
%    author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    % Extract configuration parameters
    modelName  = config.model_name;
    phaseNames = config.phaseSequence;
    x0 = config.x0;
    nx = numel(x0);
    nq = nx/2;
    nSeq = numel(phaseNames);
    traj = struct('x', [], 'u', [], 't', [], 'p', [], 'lambda', []);
    t0 = 0;
    t_max = 10;

    % Initialize ground reaction force calculation
    eventFcn = str2func([config.model_name ,config.GRFfunc]);
     
    % Initialize control parameters
    paramExtra = struct();
    paramExtra.max        = 0.2;
    paramExtra.min        = -0.2;  
    paramExtra.model_name = config.model_name;
    
    % Configure PD control parameters
    if any(strcmp(config.simControl, 'PD'))
        if strcmp(phaseNames{1}, 'Flight')
            paramExtra.alpha_target = config.simValue.alpha_target;
        else
            paramExtra.alpha_target = config.x0(nq-1);
        end

        paramExtra.k_p = config.simValue.alpha_k_p;
        paramExtra.k_d = config.simValue.alpha_k_d;
        paramExtra.tau_max = config.simValue.alpha_tau_max;
        paramExtra.idx_alpha   = nq-1;
        paramExtra.idx_dalpha  = nx-1;
    end
    
    % Configure constant control parameters
    if any(strcmp(config.simControl, 'constant')) || any(strcmp(config.simControl, 'constantAll')) 
        paramExtra.idx_dl  = numel(x0);
        paramExtra.u_alpha = config.simValue.constantInput(1);
        paramExtra.u_l     = config.simValue.constantInput(2);
    end
    
    if any(strcmp(config.simControl, 'constantHip')) 
        paramExtra.idx_dl  = numel(x0);
        paramExtra.u_alpha_const = config.simValue.constantHip;
    end

    % Configure forward simulation parameters
    if any(strcmp(config.simControl, 'forwardSimulation'))
        paramExtra.forwardSimulation = config.forwardSimulation;
        paramExtra.n_phase = 0;
        paramExtra.T = 0;
    end

    % Get model parameters and ODE settings
    p = getModelParameters(config, 0, 0);
    ode.settings = odeset('RelTol',1e-7,'AbsTol',1e-8);
    ode.solver   = @ode45;

    % Prepare control function types
    ctrlType = config.simControl;
    if ischar(ctrlType) || isstring(ctrlType)
        ctrlType = repmat({char(ctrlType)}, 1, nSeq);
    elseif iscell(ctrlType) && numel(ctrlType) == 1
        ctrlType = repmat(ctrlType, 1, nSeq);
    end

    % Simulate each phase in sequence
    for i = 1:nSeq
        phase = phaseNames{i};

        % Load phase-specific model functions
        dynamics = str2func([modelName, '.', phase, '.dynamics']);
        get_sizes = str2func([modelName, '.', phase, '.get_sizes']);
        
        % Set time step for current phase
        if numel(dt) > 1
            dt_ = dt(i);
        else
            dt_ = dt;
        end
      
        % Verify state dimensions match phase requirements
        [nx, ~, nu, ~] = get_sizes();
        assert(numel(x0) == nx, 'Initial condition x0 does not match expected size for phase: %s', phase);

        % Update forward simulation timing
        if strcmp(ctrlType(i), 'forwardSimulation')
            paramExtra.n_phase = i;
            paramExtra.T = t0;
        end

        % Create control input function
        u_fun = @(t, x) controlInput(t, x, nu, ctrlType(i), paramExtra);

        % Define phase transitions
        if nSeq > 1
            if i < nSeq
                phaseFrom = phase;
                phaseTo   = phaseNames{i+1};
            else
                phaseFrom = phase;
                phaseTo   = phaseNames{1};
            end
        elseif nSeq == 1
            phaseFrom = phase;
            phaseTo   = config.simValue.endOfSinglePhase;
        else 
            error('Invalid number of sequence phases')
        end

        % Load transition functions
        transPkg = [modelName, '.', phase, '.', phaseFrom(1), '2', phaseTo(1)];
        event_fun_raw = str2func([transPkg '.event']);
        jump_fun_raw  = str2func([transPkg '.jump']);

        % Define ODE function and event handler
        ode_fun = @(t, x) dynamics(x, u_fun(t, x), p);
        event_wrapper = @(t, x) event_general(event_fun_raw, u_fun, t, x, p);
        ode_opts = odeset('Events', event_wrapper, 'RelTol', ode.settings.RelTol, ...
                          'AbsTol', ode.settings.AbsTol);

        % Integrate dynamics
        tSpan = t0:dt_:t_max;
        if length(tSpan) > 1
            [tk, xk] = ode.solver(ode_fun, tSpan, x0, ode_opts);
        else
            tk = tSpan;
            xk = x0';
        end
        
        % Save position for specific model types during stance-to-flight transition
        if (strcmp(modelName, 'HF_noPitch_min') || strcmp(modelName, 'HF_pasHip_min') || strcmp(modelName, 'HF_compHip_min')) && strcmp(phaseTo, 'Stance')
            x_temp = xk(end,:);
        end

        % Evaluate control inputs
        uk = u_fun(tk, xk')';

        % Store trajectory data
        traj(i).x = xk(:, 1:nx);
        traj(i).u = uk;
        traj(i).t = tk - t0;
        p_opt_idx = getParameterIndex(config.optParameterNames, config);
        traj(i).p = repmat(p(p_opt_idx)', numel(tk), 1);

        % Calculate constraint forces for stance phase
        if strcmp(phase, 'Stance')
            traj(i).lambda = eventFcn(traj(i).x', traj(i).u', p)';
        else
            traj(i).lambda = zeros(size(traj(i).x, 1), 1);
        end

        % Apply jump map for phase transition
        if i < nSeq && nSeq > 1
            t0 = tk(end);
            if (strcmp(modelName, 'HF_noPitch_min') || strcmp(modelName, 'HF_pasHip_min') || strcmp(modelName, 'HF_compHip_min')) && strcmp(phase, 'Stance') && strcmp(phaseTo, 'Flight')
                x0 = jump_fun_raw(xk(end, :)', uk(end, :)', p, x_temp(1), x_temp(3), x_temp(4));
            else
                x0 = jump_fun_raw(xk(end, :)', uk(end, :)', p); 
            end
        end
    end
end

function [v, isterm, dir] = event_general(event_fun, control_fun, t, x, p)
% EVENT_GENERAL Wrapper function for phase transition events
    u = control_fun(t, x);
    v = event_fun(x, u, p);
    isterm = 1;
    dir = -1;
end

function u = controlInput(t, x, nu, ctrlType, paramExtra)
% CONTROLINPUT Generate control inputs based on specified control type
%
%   Inputs:
%     t         - time
%     x         - state vector
%     nu        - number of control inputs
%     ctrlType  - control type string
%     paramExtra- control parameters structure
%
%   Output:
%     u         - control input vector

    if strcmp(ctrlType, 'zero')
        u = zeros(nu, size(x, 2));
    elseif strcmp(ctrlType, 'PD')
        u_l = zeros(1, size(x, 2));
        u_alpha = - paramExtra.k_p * (x(paramExtra.idx_alpha, :) - paramExtra.alpha_target) - paramExtra.k_d * x(paramExtra.idx_dalpha, :);
        u_alpha = max(min(u_alpha, paramExtra.tau_max), -paramExtra.tau_max);
        u = [u_alpha; u_l];
    elseif strcmp(ctrlType, 'constant')
        u_alpha = paramExtra.u_alpha * ones(1, size(x, 2));
        u_l = zeros(1, size(x, 2));
        idx = abs(x(end, :)) < 0.7;
        u_l(idx) = paramExtra.u_l;
        u = [u_alpha; u_l];
    elseif strcmp(ctrlType, 'constantAll')
        u_alpha = paramExtra.u_alpha * ones(1, size(x, 2));
        u_l = ones(1, size(x, 2))*paramExtra.u_l;
        u = [u_alpha; u_l];
    elseif strcmp(ctrlType, 'constantHip')
        u_alpha = paramExtra.u_alpha_const * ones(1, size(x, 2));
        u_l = zeros(1, size(x, 2));
        u = [u_alpha; u_l];
    elseif strcmp(ctrlType, 'forwardSimulation')
        u_fun = paramExtra.forwardSimulation{paramExtra.n_phase};
        u = u_fun(t - paramExtra.T);
    elseif strcmp(ctrlType, 'random')
        u = paramExtra.min + (paramExtra.max - paramExtra.min) * rand(nu, size(x, 2));
    else
        error('Unknown control type: %s', ctrlType);
    end

    % Adjust output for specific model type
    if strcmp(paramExtra.model_name, 'HF_pasHip_min')
        u = u(end-nu+1:end, :);
    end
end