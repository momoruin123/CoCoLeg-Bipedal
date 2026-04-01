function traj_guess = initial_guess(config, init_stepTime, timePercentage, isParamsFlag)
modelName = config.model_name;
T = init_stepTime;
P = getModelParameters(config, 0, 0);
optParams = config.optParameterInit;
traj_guess = struct();

if strcmp(modelName, 'BF_run')
    swingFootPt_func = str2func([modelName, '.SingleStance', '.pos_SwingFoot_f']);
    nq = 5;
    nu = 4;
    N  = config.N;           % Extract segment number
    init_x = config.x0;
    L = config.operatingCond.v_avg * init_stepTime;

    % Transformation matrix
    Trans_matrix = [1 0 0 0 0;
        0 0 1 0 0;
        0 1 0 0 0;
        0 0 0 0 1;
        0 0 0 1 0;];

    % Initial guess of q
    q0_guess = init_x(1:nq);
    qM_guess = [-0.5; 0.4; 1.2; -0.5; -1];
    qT_guess = Trans_matrix * q0_guess;

    % SingleStance Phase
    N_S = N(1);
    percentage_S = timePercentage;
    T_S = percentage_S*T;
    hk = init_stepTime / N_S; % Time interval between each Segment% Initial guess of q
    t_S_init = zeros(1, N_S+1);
    % q
    q_S_init = zeros(nq, N_S+1);
    for k = 1:N_S+1
        t_k = (k-1) * hk;
        t_S_init(k) = t_k;
        q_S_init(:,k) = q0_guess + (t_k/init_stepTime)*(qM_guess - q0_guess);
    end
    % dq
    dq0_S_guess = (qM_guess-q0_guess)/T_S;
    dq_S_guess(:, :) = repmat(dq0_S_guess, 1, N_S+1);

    u_S_init = zeros(nu, N_S+1);
    traj_guess(1).x = [q_S_init; dq_S_guess]';
    traj_guess(1).u = u_S_init';
    traj_guess(1).t = t_S_init';
    traj_guess(1).Phase = 'SingleStance';

    % Flight Phase
    N_F = N(2);
    percentage_F = 1-timePercentage;
    T_F = percentage_F*T;
    hk = init_stepTime / N_F;   % Time interval between each Segment% Initial guess of q
    t_init = zeros(1, N_F+1);

    % q
    q_F_init = zeros(nq, N_F+1);
    base_init = zeros(2, N_F+1);
    stanceFootPt_T = swingFootPt_func([L; 0; q0_guess], P);
    for k = 1:N_F+1
        t_k = (k-1) * hk;
        t_init(k) = t_k;
        q_F_init(:,k) = qM_guess + (t_k/init_stepTime)*(qT_guess - qM_guess);
        base_init(:, k) = [0; 0] + (t_k/init_stepTime)*(stanceFootPt_T);
    end
    q_F_init = [base_init; q_F_init];

    % dq
    qT_guess = [stanceFootPt_T; qT_guess];
    qM_guess = [0; 0; qM_guess];
    dq0_F_guess = (qT_guess-qM_guess)/T_F;
    dq_F_guess(:,:) = repmat(dq0_F_guess, 1, N_F+1);

    u_S_init = zeros(nu, N_F+1);
    traj_guess(2).x = [q_F_init; dq_F_guess]';
    traj_guess(2).u = u_S_init';
    traj_guess(2).t = t_S_init';
    traj_guess(2).Phase = 'Flight';

    if isParamsFlag
        traj_guess(1).p = repmat(optParams, N_S+1, 1);
        traj_guess(2).p = repmat(optParams, N_F+1, 1);
    end

elseif strcmp(modelName, 'BF_min')
    nq = 5;
    nu = 4;
    N  = config.N;           % Extract segment number
    hk = init_stepTime/N;   % Time interval between each Segment% Initial guess of q
    init_x = config.x0;
    % Transformation matrix
    Trans_matrix = [1 0 0 0 0;
                    0 0 1 0 0;
                    0 1 0 0 0;
                    0 0 0 0 1;
                    0 0 0 1 0;];

    % Initial guess of q
    q0_guess = init_x(1:nq);
    qT_guess = Trans_matrix * q0_guess;

    dq0_guess = (qT_guess-q0_guess)/T;

    % linear interpolation
    t_init = zeros(1, N+1);
    q_init = zeros(nq, N+1);
    for k = 1:N+1
        t_k = (k-1)*hk;
        t_init(k) = t_k;
        q_init(:,k) = q0_guess + (t_k/init_stepTime)*(qT_guess - q0_guess);
    end

    dq_init = zeros(nq, N+1);
    dq_init(:,:) = repmat(dq0_guess, 1, N+1);

    % Initial guess of u
    u_init = zeros(nu, N+1);

    traj_guess.x = [q_init; dq_init]';
    traj_guess.u = u_init';
    traj_guess.t = t_init';
    traj_guess.Phase = 'SingleStance';

    if isParamsFlag
        traj_guess.p = repmat(optParams, N+1, 1);
    end
else
    error('Model name is incorrect or not supported');
end

end