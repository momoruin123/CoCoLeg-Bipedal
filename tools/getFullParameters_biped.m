function [fullParams, paramNames] = getFullParameters_biped(config)
% GETFULLPARAMETERS Generate complete parameter vector for hopping model
%
%   Constructs physical parameters from base values and dimensionless coefficients.
%   Supports customization through config.paramValues structure.
%
%   Inputs:
%     config - system configuration structure with optional paramValues field
%
%   Outputs:
%     fullParams - vector of physical parameter values  
%     paramNames - cell array of corresponding parameter names
%
%  author: Iskandar Khemakhem, Oussama Barhoumi, IAMS, Uni Stuttgart, 2025

    % Base dimensional parameters with default values
    m   = 32;    % Total mass [kg]
    g   = 9.81;  % Gravity [m/s²]
    l   = 0.8;   % Leg length [m]
    
    % k_h = config.paramValues.k_h;
    % k_k = config.paramValues.k_k;
    
    % Default dimensionless coefficients (fractions and scaling factors)
    coeff.m_2      = 6.8 / m;      % Thigh mass fraction
    coeff.m_3      = 3.2 / m;      % Shank mass fraction  
    coeff.m_1      = 1 - 2*coeff.m_2 - 2*coeff.m_3;  % Torso mass fraction
    coeff.theta_1  = 1.33 / (m * l^2);     % Tosro inertia scaling
    coeff.theta_2  = 0.47 / (m * l^2);     % Leg inertia scaling
    coeff.theta_3  = 0.2  / (m * l^2);      % Shank inertia scaling
    coeff.l_1      = 0.63 / l;     % Torso length
    coeff.l_2      = 0.4 / l;      % Thigh length
    coeff.d_1      = 0.315 / l;    % distance from torso CoM to hip
    coeff.d_2      = 0.2 / l;      % distance from hip to thigh CoM
    coeff.d_3      = 0.2 / l;      % distance from knee to shank CoM
    coeff.alpha0   = 0;        % Hip spring reference position
    coeff.beta0    = 0;        % Kenn spring reference position
    % coeff.xi_h     = 15.6 / (m * sqrt(g * l^3)); % Hip damping ratio
    % coeff.xi_k     = 5.48 / (m * sqrt(g * l^3));      % Knee damping ratio
    coeff.xi_h     = 2.5 / (m * sqrt(g * l^3)); % Hip damping ratio
    coeff.xi_k     = 1.1 / (m * sqrt(g * l^3));      % Knee damping ratio
    coeff.k_h      = 2 / (m * g * l);        % Hip stiffness scaling
    coeff.k_k      = 2 / (m * g * l);        % Knee stiffness scaling
    
    % Override default values if provided in config
    if isfield(config, 'paramValues') && isstruct(config.paramValues)
        pv = config.paramValues;

        % Override base dimensional parameters
        if isfield(pv, 'm'),   m   = pv.m;   end
        if isfield(pv, 'g'),   g   = pv.g;   end
        if isfield(pv, 'l'),   l   = pv.l;   end

        % Override dimensionless coefficients
        scalarFields = fieldnames(coeff);
        for i = 1:numel(scalarFields)
            fname = scalarFields{i};
            if isfield(pv, fname)
                coeff.(fname) = pv.(fname);
            end
        end
    end
    
    % Construct final parameter structure with physical units
    paramStruct.m   = m;
    paramStruct.g   = g;
    paramStruct.l   = l;

    % Convert dimensionless coefficients to physical parameters
    paramStruct.m_2      = coeff.m_2 * m;                  % Thigh mass [kg]
    paramStruct.m_3      = coeff.m_3 * m;                  % Shank mass [kg]
    paramStruct.theta_1  = coeff.theta_1 * m * l^2;        % Torso inertia [kg·m²]
    paramStruct.theta_2  = coeff.theta_2 * m * l^2;        % Thigh inertia [kg·m²]
    paramStruct.theta_3  = coeff.theta_3 * m * l^2;        % Shank inertia [kg·m²]
    paramStruct.l_1      = coeff.l_1 * l;                  % Torso length [m]
    paramStruct.l_2      = coeff.l_2 * l;                  % Thigh length [m]
    paramStruct.d_1      = coeff.d_1 * l;                  % distance from torso CoM to hip [m]
    paramStruct.d_2      = coeff.d_2 * l;                  % distance from hip to thigh CoM [m]
    paramStruct.d_3      = coeff.d_3 * l;                  % distance from knee to shank CoM [m]
    paramStruct.alpha0   = coeff.alpha0;
    paramStruct.beta0    = coeff.beta0;
    paramStruct.xi_h     = coeff.xi_h * m * sqrt(g * l^3);           % Hip damping [N·m·s/rad]
    paramStruct.xi_k     = coeff.xi_k * m * sqrt(g * l^3);           % Leg damping [N·m·s/rad]
    paramStruct.k_h      = coeff.k_h;          % Hip stiffness [N·m/rad]
    paramStruct.k_k      = coeff.k_k;          % Leg stiffness [N·m/rad]
    
    % Add any additional custom parameters not covered by standard set
    if isfield(config, 'paramValues') && isstruct(config.paramValues)
        paramFields = fieldnames(config.paramValues);
        for i = 1:numel(paramFields)
            field = paramFields{i};
            % Only add if not already in parameter structure or coefficients
            if ~isfield(paramStruct, field) && ~isfield(coeff, field)
                paramStruct.(field) = config.paramValues.(field);
            end
        end
    end

    % Return parameter vector and names for bipedal models
    if startsWith(config.model_name, 'b', 'IgnoreCase', true)
        [fullParams, paramNames] = s2v(paramStruct);
    else
        warning('This function is designed for bipedal models. Model "%s" may not be supported.', config.model_name);
        fullParams = [];
        paramNames = {};
    end
end
