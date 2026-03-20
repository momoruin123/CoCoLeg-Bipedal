function w_norm = normalize_cost(config, w, Z, P)
% Normalize cost based on total time over a hybrid sequence
%
% Inputs:
%   - config.phaseSeq: cell array of phase names (e.g., {'Flight', 'Stance', 'Flight'})
%   - config.N_array: array with number of segments per phase
%   - config.costNormalization: 'totalTime' or 'none'
%   - Z: cell array where Z{i} is the decision vector for phase i
%
% Output:
%   - w_norm: normalized cost
%
% author: Iskandar Khemakhem, Oussama Barhoumi, IAMS, Uni Stuttgart, 2025

    totalTime = 0;
    totalFlightTime = 0;
    nPhases = numel(config.phaseSequence);
    modelName = config.model_name;

    for i = 1:nPhases
        phase = config.phaseSequence{i};
        N     = config.N(i);
        Z_i   = Z{i};

        % Get size info for this phase
        get_sizes = str2func([modelName, '.', phase, '.get_sizes']);
        [nx, ~, nu, ~] = get_sizes();

        % Duration is the last element in Z{i}
        dt = Z_i(nx+nu+1);

        % Add contribution to total time
        totalTime = totalTime + N * dt;
        if strcmp(phase, 'Flight')
            totalFlightTime = totalFlightTime + N*dt;
        end
    end

    % Apply normalization
    if strcmp(config.costNormalization, 'totalTime')
        w_norm = w / totalTime;
    elseif strcmp(config.costNormalization, 'CoT_v_avg')
        w_norm = w / (config.operatingCond.v_avg*totalTime*P(1)*P(2));
        % w_norm = w/(config.operatingCond.v_avg*totalTime);
    elseif strcmp(config.costNormalization, 'FlightTime')
        w_norm = w / totalFlightTime;
    elseif strcmpi(config.costNormalization, 'none')
        w_norm = w;
    else
        error('Normalization type "%s" not available.', config.costNormalization);
    end
end
