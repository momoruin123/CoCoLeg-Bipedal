function Z = interpolate_traj2Z(config, traj, N_array)
% INTERPOLATE_TRAJ2Z Convert trajectory data to optimization decision variables
%
%   Inputs:
%     config  - system configuration structure
%     traj    - trajectory data structure array
%     N_array - number of collocation points per phase
%
%   Output:
%     Z       - vector of decision variables [states_1; controls_1; time_1; parameters_1, states_2; controls_2; time_2; parameters_2]

    % If N_array not given, try to get from config
    if nargin < 3 || isempty(N_array)
        if isfield(config, 'N')
            N_array = config.N;
        else
            error('N_array not provided and config.N is missing.');
        end
    end

    nPhases = numel(traj);
    Z = [];

    for i = 1:nPhases
        traj_i = traj(i);
        N_i    = N_array(i);

        % Duration of the phase
        T_i = traj_i.t(end);

        % Interpolate states and controls
        x_interp = interp1(traj_i.t, traj_i.x, linspace(traj_i.t(1), traj_i.t(end), N_i+1)');
        u_interp = interp1(traj_i.t, traj_i.u, linspace(traj_i.t(1), traj_i.t(end), N_i+1)');

        % Parameter vector: match order with config.optParameterNames
        if isfield(traj_i, 'p') && isfield(config, 'optParameterNames') && ~isempty(traj_i.p)
            p_interp = interp1(traj_i.t, traj_i.p, linspace(traj_i.t(1), traj_i.t(end), N_i+1)');
        else
            p_interp = [];
        end

        % Build Z_i and reshape into vector
        if config.optimizeTimeFlag
            Z_i = [x_interp, ...
                u_interp, ...
                (T_i / N_i) * ones(N_i+1, 1), ...
                p_interp]';
        else
            Z_i = [x_interp, ...
                u_interp, ...
                p_interp]';
        end
        Z_i = reshape(Z_i, [], 1);
        Z = [Z; Z_i];
    end
end