% System: hopping robot with fixed main body pitch
% phase: Stance
% Generates dynamics, constraints, and transition functions for optimization
%
% Author: Maximilian Raff, Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025
cd(fileparts(mfilename('fullpath')))
%% Floating Base Description
nx = 4;  % Number of states
nq = 2;  % Number of generalized coordinates

% Generalized coordinates and velocities
alpha  = sym('alpha','real');   % Hip angle
l      = sym('l','real');       % Leg length
dalpha = sym('dalpha','real');  % Hip angular velocity  
dl     = sym('dl','real');      % Leg extension velocity
t      = sym('t', 'real');      % Time

% Control inputs
syms u_h u_l real  % Hip torque and leg force

% State vectors
q  = [alpha l]';    % Generalized coordinates
dq = [dalpha dl]';  % Generalized velocities

%% Model Parameters
% Define symbolic physical parameters
m       = sym('m','real');      % Total mass
m_l     = sym('m_l','real');    % Leg mass
m_f     = sym('m_f','real');    % Foot mass
theta_l = sym('theta_l','real'); % Leg inertia
theta_f = sym('theta_f','real'); % Foot inertia
g       = sym('g','real');      % Gravity
l_0     = sym('l_0','real');    % Natural leg length
d_l     = sym('d_l','real');    % Leg COM offset
d_f     = sym('d_f','real');    % Foot COM offset
r_f     = sym('r_f','real');    % Foot radius
phi_0   = sym('phi_0','real');  % Neutral hip angle
xi_h    = sym('xi_h','real');   % Hip damping ratio
xi_l    = sym('xi_l','real');   % Leg damping ratio
k_h     = sym('k_h','real');    % Hip stiffness
k_l     = sym('k_l','real');    % Leg stiffness

% Gravity vector and mass distribution
g_vec = [0; -g];
m_t = m - m_f - m_l;  % Torso mass

% Damping coefficients from damping ratios
b_l = 2*xi_l*sqrt(k_l*(m_t+m_l));
b_h = 2*xi_h*sqrt(k_h*(theta_l+m_l*d_l^2+theta_f+m_f*(l_0+d_f)^2));

%% Initial state from previous phase
x0     = sym('x0','real');      % Initial x-position
alpha0 = sym('alpha0','real');  % Initial hip angle
l0     = sym('l0','real');      % Initial leg length

%% Position Reconstruction from Minimal Coordinates
% Reconstruct Cartesian positions from minimal coordinates
d0    = x0 + (phi_0 + alpha0) * r_f + l0 * sin(phi_0 + alpha0);
x_rec = d0 - r_f * (phi_0 + alpha) - l*sin(phi_0 + alpha);       % Torso x-position
y_rec = r_f + l*cos(phi_0 + alpha);                              % Torso y-position

%% Kinematics
% Compute positions and velocities of all body segments
Angle_l   = alpha + phi_0;  % Leg angle
Angle_f   = alpha + phi_0;  % Foot angle

% Center of gravity positions
CoG_t     = [x_rec; y_rec];  % Torso CoG
CoG_l     = CoG_t + d_l * [sin(Angle_l);-cos(Angle_l)];  % Leg CoG
CoG_f     = CoG_t + (l-d_f) * [sin(Angle_f);-cos(Angle_f)];  % Foot CoG

% Velocities via Jacobians
d_Angle_l = jacobian(Angle_l, q) * dq;
d_Angle_f = jacobian(Angle_f, q) * dq;
d_CoG_t   = jacobian(CoG_t, q) * dq;  % Torso velocity
d_CoG_l   = jacobian(CoG_l, q) * dq;  % Leg velocity
d_CoG_f   = jacobian(CoG_f, q) * dq;  % Foot velocity

%% Energy Formulation
% Potential energy (gravity + springs)
V_grav = CoG_t(2)*m_t*g + CoG_l(2)*m_l*g + CoG_f(2)*m_f*g;
V_spring  = 0.5*k_l*(l_0-l).^2 + 0.5*k_h*(0-alpha).^2;     % Potential Energy (due to springs)
V = V_grav + V_spring;

% Kinetic energy (translational + rotational)
T = 0.5 * (m_t     * sum(d_CoG_t.^2) + ...
           m_l     * sum(d_CoG_l.^2) + ...
           m_f     * sum(d_CoG_f.^2) + ...
           theta_l * sum(d_Angle_l.^2) + ...
           theta_f * sum(d_Angle_f.^2));

%% Euler-Lagrange Equations
% M * ddq + C * dq + G = 0
% Derive mass matrix (M), Coriolis matrix (C), and gravity vector (G) using
% Euler-Lagrange formalism
% only V_grav since the springs are dealt with via the contact matricies
[M, C, G] = eulerLagrange(T, V_grav, q, dq); 

% Spring and damping forces
Fspring_l = k_l*(l_0-l)+b_l*(0-dl); 
Fspring_h = k_h*(0-alpha)+b_h*(0-dalpha); 

% Input mapping matrices
B =  [[1;0],...
     [0;1]];

% Passive spring torques
tau_J = B*[Fspring_h; Fspring_l];

%% Continuous Dynamics
% Compute acceleration: ddq = M^(-1) * (tau - C*dq - G)
acc = M \ (-C-G+tau_J+B*[u_h;u_l]);
% State derivative
f_S = [dq; acc(1:nq)]; 

%% Event and Jump Functions for Phase Transitions
% S2F: Lift-off event (transition from stance to flight): net vertical force becomes zero
jac_CoG_f = [jacobian(d_CoG_f, q), jacobian(d_CoG_f, dq)];
jac_CoG_t = [jacobian(d_CoG_t, q), jacobian(d_CoG_t, dq)];
jac_CoG_l = [jacobian(d_CoG_l, q), jacobian(d_CoG_l, dq)];
% Vertical acceleration Jacobian for all bodies
threeBodiesAccjac_Y = [jac_CoG_t(2,:); jac_CoG_l(2,:); jac_CoG_f(2,:)];
dd_y = threeBodiesAccjac_Y * f_S;

event_lo = [m_t, m_l, m_f] * (dd_y + g);
event_lo = simplify(event_lo);

% Lift-off jump map (full and minimal versions)
jump_lo  = [x_rec ; y_rec ; q ; d_CoG_t ;dq];
% define an additional jump map that ignores the x direction to use in
% casadi context. When optimizing in a flight stance cylce, we don't care
% about the final x position in the jump function as this is only needed in
% the periodicity constraints which doesn't look at the x position
jump_noX  = [0; y_rec; q; d_CoG_t; dq];

% S2S: Nadir event (minimum leg length)
event_nadir = -dl;    
jump_nadir  = [q; dq];

%% Compute positive mechanical work
% Smooth approximation of mechanical power
posMechWork = smoothRampFcn(dl*u_l + dalpha*u_h, 0.3);

%% Total Energy
E = T + V;

%% parameters vector 
% parameters to be optimized should always placed AT THE END OF THE VECTOR
p = [m; g; l_0; m_l; m_f; theta_l; theta_f; d_l; d_f; r_f; phi_0; xi_h; xi_l; k_h; k_l];  % k_l is the only param to be optimized on here

%% Automatic Function Generation
% Generate MATLAB functions for the dynamics and their derivatives
x = [q; dq];   % State vector
u = [u_h;u_l]; % input vector

% Generate MATLAB functions and export them to files
matlabFunction(f_S, 'File', 'dynamics', 'Vars', {x, u, p});
matlabFunction(M, 'File', 'massMatrix', 'Vars', {x, p});
matlabFunction(-C-G+tau_J+B*u, 'File', 'forces', 'Vars', {x, u, p});
matlabFunction(E, 'File', 'energyLevel', 'Vars', {x, u, p});
matlabFunction(posMechWork, 'File', 'posMechWork', 'Vars', {x, u});

%% Create packages for transitions
% Define package folder names (in current folder)
packageNames = {'S2F', 'S2S'};

% Create package folders in current directory if they don't exist
for i = 1:length(packageNames)
    packageName = packageNames{i};
    packageFolder = ['+', packageName];

    if ~exist(packageFolder, 'dir')
        mkdir(packageFolder);
    end
end

% Generate the MATLAB functions in the appropriate package folders
matlabFunction(event_lo, 'File', fullfile(['+', packageNames{1}], 'event'), 'Vars', {x, u, p});
matlabFunction(jump_lo, 'File', fullfile(['+', packageNames{1}], 'jump'), 'Vars', {x, u, p, x0, alpha0, l0});
matlabFunction(jump_noX, 'File', fullfile(['+', packageNames{1}], 'jump_noX'), 'Vars', {x, u, p});
matlabFunction(threeBodiesAccjac_Y, 'File',  fullfile(['+', packageNames{1}], 'threeBodiesAccjac_Y'), 'Vars', {x, p});
matlabFunction([m_t, m_l, m_f], 'File', fullfile(['+', packageNames{1}], 'masses'), 'Vars', {p});
matlabFunction(event_nadir, 'File', fullfile(['+', packageNames{2}], 'event'), 'Vars', {x, u, p});
matlabFunction(jump_nadir, 'File', fullfile(['+', packageNames{2}], 'jump'), 'Vars', {x, u, p});

%% Generate Model Dimension Functions
% Define sizes in a function
nx = numel(x);
nu = numel(u);
np = numel(p);
nq = numel(q);

% Create a function file that returns these sizes
fid = fopen('get_sizes.m', 'w');
fprintf(fid, 'function [nx, nq, nu, np] = get_sizes()\n');
fprintf(fid, 'nx = %d;\n', nx);
fprintf(fid, 'nq = %d;\n', nq);
fprintf(fid, 'nu = %d;\n', nu);
fprintf(fid, 'np = %d;\n', np);
fprintf(fid, 'end\n');
fclose(fid);

%% Create a function to define the needed parameters
% This function is only needed once and is valid for the entire model
% Extract names.
paramNamesAll = arrayfun(@char, p, 'UniformOutput', false);

% === Create function in parent folder ===
fid = fopen('../getParameterNames.m', 'w');

% Write function header
fprintf(fid, 'function paramNames = getParameterNames(config)\n');
fprintf(fid, '%%GETPARAMETERNAMES_HOPPER Returns the parameter names for the Hopper model\n\n');

% Write initial full parameter list
fprintf(fid, 'paramNamesAll = {\n');
for i = 1:length(paramNamesAll)
    fprintf(fid, '    ''%s'';\n', paramNamesAll{i});
end
fprintf(fid, '};\n\n');

% Start handling the config.optParameterNames
fprintf(fid, '%% Reorder if config.optParameterNames is given\n');
fprintf(fid, 'if isfield(config, ''optParameterNames'') && ~isempty(config.optParameterNames)\n');
fprintf(fid, '    optNames = config.optParameterNames;\n');
fprintf(fid, '    %% Keep only those optimization names that exist\n');
fprintf(fid, '    optNames = optNames(ismember(optNames, paramNamesAll));\n');
fprintf(fid, '    %% Remove optimization names from main list\n');
fprintf(fid, '    paramNames = setdiff(paramNamesAll, optNames, ''stable'');\n');
fprintf(fid, '    %% Append optimization names at the end\n');
fprintf(fid, '    paramNames = [paramNames; optNames(:)];\n');
fprintf(fid, 'else\n');
fprintf(fid, '    paramNames = paramNamesAll;\n');
fprintf(fid, 'end\n\n');

% End of function
fprintf(fid, 'end\n');

% Close file
fclose(fid);

%% Define the function that indicates the position of the minimal states and inputs in the full states
idx_x = [4;5;9;10];  
idx_u = [1;2];  

% === Create function file in parent folder ===
fid = fopen('get_classIndices.m', 'w');  % Save in parent folder

% Write function header
fprintf(fid, 'function [idx_x, idx_u] = get_classIndices()\n');
fprintf(fid, '%%GET_CLASSINDICES Returns the indices of the optimal states and inputs in the full states\n\n');

% Write state indices
fprintf(fid, '%% State indices\n');
fprintf(fid, 'idx_x = [');
for i = 1:length(idx_x)-1
    fprintf(fid, '%d; ', idx_x(i));
end
fprintf(fid, '%d', idx_x(end));
fprintf(fid, '];\n\n');

% Write input indices
fprintf(fid, '%% Input indices\n');
fprintf(fid, 'idx_u = [');
for i = 1:length(idx_u)-1
    fprintf(fid, '%d; ', idx_u(i));
end
fprintf(fid, '%d', idx_u(end));
fprintf(fid, '];\n\n');

% End function
fprintf(fid, 'end\n');

% Close file
fclose(fid);

%% Helper Functions
function [M, C, G] = eulerLagrange(T, V, q, dq)
    % Compute the mass matrix (M), Coriolis forces (C), and gravity vector (G)
    % using the Euler-Lagrange formalism.

    % Partial derivatives
    dT_dq   = jacobian(T, q).';  % Partial derivative of kinetic energy w.r.t. q
    dV_dq   = jacobian(V, q).';  % Partial derivative of potential energy w.r.t. q
    dT_dqdt = jacobian(T, dq).'; % Partial derivative of Lagrangian w.r.t. dq

    % Time derivative of dL/dqdt
    dd_T_dqdt2  = jacobian(dT_dqdt, dq);
    d_dTdqdt_dq = jacobian(dT_dqdt, q);

    % Assign matrices
    M = simplify(dd_T_dqdt2);  % Mass matrix
    C = -simplify(expand(dT_dq - d_dTdqdt_dq * dq));  % Coriolis vector
    G = simplify(expand(dV_dq));  % Gravity vector
end

function [W, W_dot] = computeContactMatrix(con, q, dq)
    % Compute the contact projection matrix W and its time derivative W_dot
    % for a given constraint 'con'.
    
    W = jacobian(con, q)';  % Contact matrix
    W_dot = W;

    % Compute time derivative of each column of W
    for j = 1:size(W, 2)
        W_dot(:, j) = jacobian(W(:, j), q) * dq;  % Time derivative of contact matrix
    end
end
