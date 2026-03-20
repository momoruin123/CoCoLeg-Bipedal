function traj = interpolate_Z2traj(config, Z, N_array, interpolateFlag, dt, dxFlag)
% INTERPOLATE_Z2TRAJ Convert flattened decision variables to trajectory with optional interpolation
%
%   Converts flattened decision variables from optimization into structured
%   trajectory data, with optional Hermite-Simpson interpolation for finer
%   resolution and derivative computation.
%
%   Inputs:
%     config         - system configuration structure
%     Z              - flattened decision variable vector
%     N_array        - number of segments per phase
%     interpolateFlag- enable Hermite-Simpson interpolation (default: false)
%     dt             - time step for interpolation (required if interpolateFlag=true)
%     dxFlag         - compute state derivatives (only with interpolation)
%
%   Output:
%     traj           - trajectory structure with fields:   
%                      x: states, u: controls, t: time, p: parameters
%                      lambda: ground reaction forces, dx: state derivatives
%
%   author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025


    % Initialize ground reaction force calculation function
    eventFcn = str2func([config.model_name ,config.GRFfunc]);
    
    % Handle optional inputs with defaults
    if nargin < 3 || isempty(N_array)
        if isfield(config, 'N')
            N_array = config.N;
        else
            error('N_array not provided and config.N is missing.');
        end
    end

    if nargin < 4
        interpolateFlag = false;
    end

    if nargin == 4
        error('Interpolation time step must be specified when interpolateFlag is true');
    end

    if nargin < 6
        dxFlag = false;
    end

    if dxFlag && ~interpolateFlag
        warning('State derivative computation only available with interpolation enabled');
    end

    % Initialize trajectory structure and counters
    nPhases = numel(config.phaseSequence);
    traj = struct('x', {}, 'u', {}, 't', {}, 'p', {} , 'lambda' , {}, 'Phase', {});
    na = numel(config.optParameterNames) + config.optimizeTimeFlag;
    idx_Z = 1;

    % Process each phase
    for i = 1:nPhases
        N_i = N_array(i);
        phase_name = config.phaseSequence{i};
        
        traj(i).Phase = phase_name;
        
        % Load phase-specific model functions
        get_sizes = str2func([config.model_name, '.', phase_name, '.get_sizes']);
        dynamics  = str2func([config.model_name, '.', phase_name, '.dynamics']);
        [nx, ~, nu, ~] = get_sizes();
        
        % Extract decision variables for current phase
        rows_per_point = nx + nu + na;
        total_rows = rows_per_point * (N_i + 1);
        Z_i = Z(idx_Z : idx_Z + total_rows - 1);
        Z_i = reshape(Z_i, [], N_i+1);
        
        % Parse variables from decision vector
        X = Z_i(1:nx, :)';                    % States
        U = Z_i(nx+1:nx+nu, :)';              % Controls  
        h = Z_i(nx+nu+1, 1);                  % Time step
        p = Z_i(end-na+2:end, 1)';            % Optimized parameters
        
        % Construct full parameter vector
        parameters = getModelParameters(config, 0, 1);
        p_full = [parameters; p'];
        
        % Create coarse time grid
        gridTime = linspace(0, N_i*h, N_i+1)';

        if interpolateFlag
            % Create fine time grid for interpolation
            t_array = (0:dt:N_i*h)';
            if t_array(end) < N_i*h
                t_array = [t_array; N_i*h];  % Ensure endpoint inclusion
            end
            
            if dt < h
                % Hermite-Simpson interpolation for fine resolution
                closest_indices = zeros(length(t_array),1);
                closest_indices(1) = 1;
                
                % Find closest coarse grid indices for each fine grid point
                for j = 2:numel(t_array)
                    indices_below_t = find(gridTime < t_array(j));
                    closest_indices(j) = indices_below_t(end);
                end
                
                % Compute normalized time within each interval
                tau = t_array(1:length(closest_indices)) - gridTime(closest_indices);
    
                % Extract variables at interval boundaries
                u_k   = U(closest_indices,:);
                u_k1  = U(closest_indices+1,:);
                x_k   = X(closest_indices,:);
                x_k1  = X(closest_indices+1,:);
                
                % Compute dynamics at boundaries
                f_k   = dynamics(x_k', u_k', p_full)';
                f_k1  = dynamics(x_k1', u_k1', p_full)';

                % Compute midpoint values using specified schemes
                [x_mid, u_mid, f_mid] = computeMidpointValues(config, x_k, x_k1, u_k, u_k1, f_k, f_k1, h, dynamics, p_full);
                
                 % Interpolate using specified schemes
                [traj(i).x, traj(i).u, traj_dx] = interpolateWithScheme(config, x_k, x_k1, u_k, u_k1, f_k, f_mid, f_k1, tau, h, dxFlag);

                % incude auxialliary variables
                traj(i).t = t_array;
                traj(i).p = repmat(p, numel(t_array), 1);
                
                % Compute state derivatives if requested
                if dxFlag
                    traj(i).dx = traj_dx;
                end
     
                % Compute ground reaction forces for stance phase
                if strcmp(phase_name, 'Stance')
                    traj(i).lambda = eventFcn(traj(i).x', traj(i).u', p_full)';
                else
                    traj(i).lambda = zeros(size(traj(i).x, 1), 1);
                end

            else
                % Coarse grid sampling with linear interpolation
                traj_x = zeros(length(t_array), nx);
                traj_u = zeros(length(t_array), nu);
                
                for j = 1:length(t_array)
                    t = t_array(j);
                    % Find bounding indices in coarse grid
                    idx_left = find(gridTime <= t, 1, 'last');
                    if idx_left == N_i+1
                        idx_left = N_i;  % Prevent index overflow
                    end
                    idx_right = idx_left + 1;
            
                    % Linear interpolation weights
                    tL = gridTime(idx_left);
                    tR = gridTime(idx_right);
                    alpha = (t - tL) / (tR - tL + eps);  % eps prevents division by zero
                    
                    % Interpolate states and controls
                    traj_x(j,:) = (1-alpha)*X(idx_left,:) + alpha*X(idx_right,:);
                    traj_u(j,:) = (1-alpha)*U(idx_left,:) + alpha*U(idx_right,:);
                end
            
                traj(i).x = traj_x;
                traj(i).u = traj_u;
                traj(i).t = t_array;
                traj(i).p = repmat(p, numel(t_array), 1);
            
                % Compute ground reaction forces
                if strcmp(phase_name, 'Stance')
                    traj(i).lambda = eventFcn(traj_x', traj_u', p_full)';
                else
                    traj(i).lambda = zeros(size(traj_x, 1), 1);
                end
            end

        else
            % No interpolation - use coarse grid directly
            traj(i).x = X;
            traj(i).u = U;
            traj(i).t = gridTime;
            traj(i).p = repmat(p, numel(gridTime), 1);

            % Compute ground reaction forces
            if strcmp(phase_name, 'Stance')
                traj(i).lambda = eventFcn(traj(i).x', traj(i).u', p_full)';
            elseif strcmp(phase_name, 'SingleStance')
                lambda = [];
                for j = 1:size(X,1)
                    lambda_j = BF_min.SingleStance.compute_lambda(X(j,:)', U(j,:)', p_full);
                    lambda = [lambda;lambda_j'];
                end
                % [~, lambda] = BF_min.SingleStance.S2S.jump(traj(i).x', traj(i).u', p_full);
                traj(i).lambda = lambda;
            else
                traj(i).lambda = zeros(size(traj(i).x, 1), 1);
            end
        end
        
        % Advance index in decision variable vector
        idx_Z = idx_Z + total_rows;
    end
end

function [x_mid, u_mid, f_mid] = computeMidpointValues(config, x_k, x_k1, u_k, u_k1, f_k, f_k1, h, dynamics, p_full)
% COMPUTEMIDPOINTVALUES Compute midpoint values using specified collocation and interpolation schemes
%
%   Inputs:
%     config     - system configuration structure
%     x_k, x_k1 - states at interval boundaries
%     u_k, u_k1 - controls at interval boundaries
%     f_k, f_k1 - dynamics at interval boundaries
%     h          - time step
%     dynamics   - dynamics function handle
%     p_full     - full parameter vector
%
%   Outputs:
%     x_mid      - midpoint states
%     u_mid      - midpoint controls
%     f_mid      - midpoint dynamics

    % Extract schemes from config with defaults
    if isfield(config, 'collocationScheme')
        collocationScheme = config.collocationScheme;
    else
        collocationScheme = 'hermiteSimpson';
    end

    if isfield(config, 'inputInterpolation')
        inputInterpolation = config.inputInterpolation;
    else
        inputInterpolation = 'piecewiseLinear';
    end

    % Compute midpoint controls based on interpolation scheme
    switch inputInterpolation
        case 'piecewiseLinear'
            u_mid = (u_k + u_k1)/2;
        case 'piecewiseConstant'
            u_mid = u_k;  % Use left endpoint value
        otherwise
            error('Unknown input interpolation scheme: %s', config.inputInterpolation);
    end

    % Compute midpoint states based on collocation scheme
    switch collocationScheme
        case 'hermiteSimpson'
            % Hermite-Simpson midpoint: x_mid = 0.5*(x_k+x_k1) + (h/8)*(f_k-f_k1)
            x_mid = 0.5*(x_k + x_k1) + (h/8)*(f_k - f_k1);
        case 'trapezoidal'
            % Trapezoidal midpoint: simple average
            x_mid = 0.5*(x_k + x_k1);
        otherwise
            error('Unknown collocation scheme: %s', config.collocationScheme);
    end

    % Compute dynamics at midpoint
    f_mid = dynamics(x_mid', u_mid', p_full)';
end

function [traj_x, traj_u, traj_dx] = interpolateWithScheme(config, x_k, x_k1, u_k, u_k1, f_k, f_mid, f_k1, tau, h, dxFlag)
% INTERPOLATEWITHSCHEME Interpolate states and controls using specified schemes
%
%   Inputs:
%     config  - system configuration structure
%     x_k, x_k1 - states at interval boundaries
%     u_k, u_k1 - controls at interval boundaries
%     f_k, f_mid, f_k1 - dynamics at boundaries and midpoint
%     tau     - normalized time within interval [0, h]
%     h       - time step
%     dxFlag  - flag to compute state derivatives
%
%   Outputs:
%     traj_x  - interpolated states
%     traj_u  - interpolated controls
%     traj_dx - interpolated state derivatives (if dxFlag=true)

    % Extract schemes from config with defaults
    if isfield(config, 'collocationScheme')
        collocationScheme = config.collocationScheme;
    else
        collocationScheme = 'hermiteSimpson';
    end

    if isfield(config, 'inputInterpolation')
        inputInterpolation = config.inputInterpolation;
    else
        inputInterpolation = 'piecewiseLinear';
    end

    % Interpolate controls based on scheme
    switch inputInterpolation
        case 'piecewiseLinear'
            traj_u = u_k + (u_k1 - u_k) .* (tau / h);
        case 'piecewiseConstant'
            traj_u = u_k;  % Constant within each interval
        otherwise
            error('Unknown input interpolation scheme: %s', config.inputInterpolation);
    end

    % Interpolate states based on collocation scheme
    switch collocationScheme
        case 'hermiteSimpson'
            % Cubic Hermite interpolation for Hermite-Simpson
            traj_x = x_k + f_k.*tau + 0.5*(-3*f_k + 4*f_mid-f_k1).*(tau.^2)/h + (1/3)*(2*f_k-4*f_mid+2*f_k1).*(tau.^3)/(h^2);
            
            if dxFlag
                % Derivative of cubic interpolation
                traj_dx = f_k + (-3*f_k+4*f_mid-f_k1).*(tau/h) + (2*f_k-4*f_mid + 2*f_k1).*((tau/h).^2);
            else
                traj_dx = [];
            end

        case 'trapezoidal'
            % Linear interpolation for trapezoidal method
            traj_x = x_k + (x_k1 - x_k) .* (tau / h);
            
            if dxFlag
                % Constant derivative for linear interpolation
                traj_dx = (x_k1 - x_k) / h;
            else
                traj_dx = [];
            end

        otherwise
            error('Unknown collocation scheme: %s', config.collocationScheme);
    end
end