function [fullParams, paramNames] = getFullParameters_hopper(config)
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
    m   = 1;    % Total mass [kg]
    g   = 1;    % Gravity [m/s²]
    l_0 = 1;    % Rest leg length [m]
    
    % Default dimensionless coefficients (fractions and scaling factors)
    coeff.m_f      = 0.05;     % Foot mass fraction
    coeff.m_l      = 0.1;      % Leg mass fraction  
    coeff.m_t      = m - coeff.m_f - coeff.m_l;  % Torso mass fraction
    coeff.theta_f  = 0.002*2;  % Foot inertia scaling
    coeff.theta_l  = 0.002*2;  % Leg inertia scaling
    coeff.theta_t  = 0.4;      % Torso inertia scaling
    coeff.k_h      = 5;        % Hip stiffness scaling
    coeff.xi_l     = 0.2 * sqrt(2); % Leg damping ratio
    coeff.xi_h     = 0.2;      % Hip damping ratio
    coeff.r_f      = 0.05;     % Foot radius fraction
    coeff.d_l      = 0.25;     % Leg COM fraction
    coeff.d_f      = 0.25;     % Foot COM fraction
    coeff.phi_0    = 0;        % Neutral hip angle [rad]
    coeff.k_l      = 20;       % Leg stiffness scaling
    
    % Override default values if provided in config
    if isfield(config, 'paramValues') && isstruct(config.paramValues)
        pv = config.paramValues;
    
        % Override base dimensional parameters
        if isfield(pv, 'm'),   m   = pv.m;   end
        if isfield(pv, 'g'),   g   = pv.g;   end
        if isfield(pv, 'l_0'), l_0 = pv.l_0; end
    
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
    paramStruct.l_0 = l_0;
    
    % Convert dimensionless coefficients to physical parameters
    paramStruct.m_f      = coeff.m_f * m;                    % Foot mass [kg]
    paramStruct.m_l      = coeff.m_l * m;                    % Leg mass [kg]
    paramStruct.m_t      = coeff.m_t * m;                    % Torso mass [kg]
    paramStruct.theta_f  = coeff.theta_f * m * l_0^2;        % Foot inertia [kg·m²]
    paramStruct.theta_l  = coeff.theta_l * m * l_0^2;        % Leg inertia [kg·m²]
    paramStruct.theta_t  = coeff.theta_t * m * l_0^2;        % Torso inertia [kg·m²]
    paramStruct.k_h      = coeff.k_h * m * g * l_0;          % Hip stiffness [N·m/rad]
    paramStruct.r_f      = coeff.r_f * l_0;                  % Foot radius [m]
    paramStruct.d_l      = coeff.d_l * l_0;                  % Leg COM position [m]
    paramStruct.d_f      = coeff.d_f * l_0;                  % Foot COM position [m]
    paramStruct.phi_0    = coeff.phi_0;                      % Neutral hip angle [rad]
    paramStruct.xi_h     = coeff.xi_h * m * g * l_0;         % Hip damping [N·m·s/rad]
    paramStruct.xi_l     = coeff.xi_l * m * sqrt(g / l_0);   % Leg damping [N·s/m]
    paramStruct.k_h      = coeff.k_h * m * g * l_0;          % Hip stiffness [N·m/rad]
    paramStruct.k_l      = coeff.k_l * m * g / l_0;          % Leg stiffness [N/m]

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

    % Return parameter vector and names for hopping models
    if startsWith(config.model_name, 'h', 'IgnoreCase', true)
        [fullParams, paramNames] = s2v(paramStruct);
    else
        warning('This function is designed for hopping models. Model "%s" may not be supported.', config.model_name);
        fullParams = [];
        paramNames = {};
    end
end
