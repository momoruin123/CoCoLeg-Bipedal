function [r, constraintName] = init_operating_residuals(config, Z)
%INIT_OPERATING_RESIDUALS Creates residuals for initial conditions and operating point
%    author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    import casadi.*
    [nx, nq, nu, np] = HP_min.Flight.get_sizes(); 

    % Extract first phase if Z is cell array
    if iscell(Z)
        Z = Z{1};
    end

    x_init = Z(1:nx);

    % Initialize residuals array
    r = MX.zeros(2*nx);
    constraintName = cell(2*nx,1);
    idx_r = 0;

    % Hopping height constraint: initial height = operating height
    r(idx_r+1:idx_r+numel(config.operatingCond.idxHP)) = x_init(config.operatingCond.idxHP) - config.operatingCond.HP;
    constraintName(idx_r+1:idx_r+numel(config.operatingCond.idxHP)) = {'initialization and operating point'};
    idx_r = idx_r + numel(config.operatingCond.idxHP);

    % Truncate to actual size
    r = r(1:idx_r);
    constraintName = constraintName(1:idx_r);
end