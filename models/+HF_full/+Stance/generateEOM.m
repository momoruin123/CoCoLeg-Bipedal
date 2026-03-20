% System: hopping robot with 5 DoF
% phase: Stance
% Generates dynamics, constraints, and transition functions for optimization
%
% Author: Maximilian Raff, Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025
cd(fileparts(mfilename('fullpath')))
%% Floating Base description
nx = 10;  % Number of states
nq = 5;

% Generalized coordinates
x      = sym('x','real');
y      = sym('y','real');
alpha  = sym('alpha','real');
phi    = sym('phi','real');
l      = sym('l','real');
% Generalized velocities
dx     = sym('dx','real');
dy     = sym('dy','real');
dalpha = sym('dalpha','real');
dphi   = sym('dphi','real');
dl     = sym('dl','real');

t      = sym('t', 'real');

syms u_h u_l real  % Input variable (control)

q  = [x y phi alpha l]'; % Generalized coordinates vector
dq = [dx dy dphi dalpha dl]'; % Generalized velocities vector

%% Model Parameters
% Define symbolic variables for the model's physical parameters
m       = sym('m','real');    % mass of upper body / torso
m_l     = sym('m_l','real');    % mass of leg / lower body
m_f     = sym('m_f','real');    % mass of foot / lower body
theta_t = sym('theta_t','real');
theta_l = sym('theta_l','real');
theta_f = sym('theta_f','real');
g       = sym('g','real');      % gravity
l_0     = sym('l_0','real');    % natural length of spring
d_l     = sym('d_l','real');
d_f     = sym('d_f','real');
r_f     = sym('r_f','real');       % radius of foot
k_l     = sym('k_l','real');       % leg spring stiffness
k_h     = sym('k_h','real');       % hip spring stiffness
xi_l    = sym('xi_l','real');      % leg spring damping
xi_h    = sym('xi_h','real');      % hip spring damping

% Gravitational force vector in terms of incline gamma
g_vec = [0; -g];

% Define damping ratios
m_t = m - m_f - m_l;
b_l = 2*xi_l*sqrt(k_l*(m_t+m_l));
b_h = 2*xi_h*sqrt(k_h*(theta_l+m_l*d_l^2+theta_f+m_f*(l_0+d_f)^2));

%% constraints
%  constraint functionals
syms x_c alpha_c phi_c l_c real  % initial contact point
con_S = [x+r_f*(phi+alpha)+l*sin(phi+alpha)-x_c-(alpha_c+phi_c)*r_f-l_c*sin(alpha_c+phi_c);...
         y-l*cos(phi+alpha)-r_f];   % constraint during stance
% Compute contact matrices (W) and their time derivatives (W_dot)
[W_S, W_S_dot] = computeContactMatrix(con_S, q, dq);  % rolling contact

%% Kinematics
% Centers of Gravity (CoG) and velocities computed via Jacobians
% Positions of the centers of gravity (CoG)
Angle_t   = phi;
Angle_l   = Angle_t+alpha;
Angle_f   = Angle_t+alpha;
CoG_t     = [x;y];
CoG_l     = CoG_t + d_l * [sin(Angle_l);-cos(Angle_l)];
CoG_f     = CoG_t + (l-d_f) * [sin(Angle_f);-cos(Angle_f)];
% Velocities of the centers of gravity
d_Angle_t = jacobian(Angle_t, q) * dq;
d_Angle_l = jacobian(Angle_l, q) * dq;
d_Angle_f = jacobian(Angle_f, q) * dq;
d_CoG_t   = jacobian(CoG_t, q) * dq;  % Velocity of the torso
d_CoG_l   = jacobian(CoG_l, q) * dq;  % Velocity of the leg
d_CoG_f   = jacobian(CoG_f, q) * dq;  % Velocity of the foot

%% Energies
% Define potential and kinetic energy expressions

% Potential Energy (due to gravity)
V_grav = CoG_t(2)*m_t*g + CoG_l(2)*m_l*g + CoG_f(2)*m_f*g;
% Potential Energy (due to springs)
V_spring  = 0.5*k_l*(l_0-l).^2 + 0.5*k_h*(0-alpha).^2;
% Kinetic Energy:         
T = 0.5 * (m_t     * sum(d_CoG_t.^2) + ...
           m_l     * sum(d_CoG_l.^2) + ...
           m_f     * sum(d_CoG_f.^2) + ...
           theta_t * sum(d_Angle_t.^2) + ...
           theta_l * sum(d_Angle_l.^2) + ...
           theta_f * sum(d_Angle_f.^2));
V = V_grav + V_spring;

%% Euler-Lagrange Equations
% M * ddq + C * dq + G = 0
% Derive mass matrix (M), Coriolis matrix (C), and gravity vector (G) using
% Euler-Lagrange formalism
% only V_grav since the springs are dealt with via the contact matricies
[M, C, G] = eulerLagrange(T, V_grav, q, dq); 

% Spring Forces
Fspring_l = k_l*(l_0-l)+b_l*(0-dl); 
Fspring_h = k_h*(0-alpha)+b_h*(0-dalpha); 

B = [[0;0;0;1;0],...
     [0;0;0;0;1]];

tau_J = B*[Fspring_h; Fspring_l];

%% Compute Continuos Dynamics
acc = [M -W_S; W_S' zeros(2)] \ [-C-G+tau_J+B*[u_h;u_l];-W_S_dot'*dq];
lambda_S = acc(nq+1:end);
f_S = [dq; acc(1:nq)];

%% Define event and mapping lift-off
event_lo = lambda_S(2); % transition into flight
plus = [M -W_S; W_S' zeros(2)] \ [M*dq;0;0]; 
jump_lo  = [q; plus(1:5)];

%% Define event and mapping nadir
event_nadir = -dl;    
jump_nadir  = [q; dq];

%% Compute positive mechanical work
posMechWork = smoothRampFcn(dl*u_l + dalpha*u_h, 0.3);

%% Energy level
E = T + V;
%% parameters vector
p = [m; g; l_0; m_l; m_f; theta_t;theta_l; theta_f; d_l; d_f; r_f; xi_h; xi_l; k_h; k_l];
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
matlabFunction(W_S, 'File', 'contactJacobian', 'Vars', {x, p});
matlabFunction(W_S_dot, 'File', 'contactJacobianDerivative', 'Vars', {x, p});

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
matlabFunction(jump_lo, 'File', fullfile(['+', packageNames{1}], 'jump'), 'Vars', {x, u, p});
matlabFunction(event_nadir, 'File', fullfile(['+', packageNames{2}], 'event'), 'Vars', {x, u, p});
matlabFunction(jump_nadir, 'File', fullfile(['+', packageNames{2}], 'jump'), 'Vars', {x, u, p});

%% Define sizes in a function
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
% Extract names
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

%% Define the function that indicates the position of the optimal states and inputs in the full states
idx_x = [1;2;3;4;5;6;7;8;9;10];  
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
