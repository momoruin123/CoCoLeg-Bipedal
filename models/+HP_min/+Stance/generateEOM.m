% System: Harmonic oscillator, which corresponds to a hopper in hopping in
% place and in Stance phase
cd(fileparts(mfilename('fullpath')))
%% Floating Base description
nx = 2;  % Number of states
nq = 1;
% Generalized coordinates
l      = sym('l','real');
% Generalized velocities
dl     = sym('dl','real');
% time
t      = sym('t', 'real');

syms u real  % Input variable (control)

q  = [l]'; % Generalized coordinates vector
dq = [dl]'; % Generalized velocities vector

%% Model Parameters
% Define symbolic variables for the model's physical parameters
m       = sym('m','real');    % mass of upper body / torso
m_l     = sym('m_l','real');    % mass of leg / lower body
m_f     = sym('m_f','real');    % mass of foot / lower body
g       = sym('g','real');      % gravity
l_0     = sym('l_0','real');    % natural length of spring
d_l     = sym('d_l','real');
d_f     = sym('d_f','real');
r_f     = sym('r_f','real');      % radius of foot
k_l     = sym('k_l','real');      % spring stiffness
xi_l    = sym('xi_l','real');      % spring damping

% Gravitational force vector in terms of incline gamma
g_vec = [0; -g];

% Define damping ratios
m_t = m - m_f - m_l;
m_tl = m_t + m_l;
b_l = 2*xi_l*sqrt(k_l*(m_t+m_l));


%% Kinematics
% Centers of Gravity (CoG) and velocities computed via Jacobians
% Positions of the centers of gravity (CoG)
CoG_t = r_f + l;
CoG_l = CoG_t - d_l;

CoG_tl = (m_t * CoG_t + m_l * CoG_l) / m_tl;
% Velocities of the centers of gravity
d_CoG_tl   = jacobian(CoG_tl, q) * dq;  % Velocity of the torso

%% Energies
% Define potential and kinetic energy expressions

% Potential Energy (due to gravity)
V_grav = CoG_tl*m_tl*g ;
% Potential Energy (due to springs)
V_spring  = 0.5*k_l*(l_0-l).^2; % + 0.5*k_h*(0-alpha).^2;
% Kinetic Energy:         
T = 0.5 * m_tl  * sum(d_CoG_tl.^2);
V = V_grav + V_spring;

%% Euler-Lagrange Equations
% M * ddq + C * dq + G = 0
% Derive mass matrix (M), Coriolis matrix (C), and gravity vector (G) using
% Euler-Lagrange formalism
% only V_grav since the springs are dealt with via the contact matricies
[M, C, G] = eulerLagrange(T, V_grav, q, dq); 

% Spring Forces
Fspring_l = k_l*(l_0-l)+b_l*(0-dl); 

%% Contact Projection Matrices
% Compute Continuos Dynamics
acc = M\(-C-G+Fspring_l+u);
f = [dq;acc(1:nq)];

%% organize in A, B and D 
x = [q; dq];   % State vector
A_lin = jacobian(f, x);
B_lin = jacobian(f, u);
D_lin = f - A_lin*x - B_lin*u;
D_lin = simplify(D_lin);

%% Define event and mapping lift-off
event_lo = Fspring_l + B*u + m_f*g;
jump_lo  = [l+r_f; l; dl; dl];

%% Define event and mapping nadir
event_nadir = -dl;    
jump_nadir  = [q; dq];

%% Compute positive mechanical work
posMechWork = smoothRampFcn(dl*u, 0.3);

%% Energy level
E = T + V;
%% parameters vector
% p = [m; g; l_0; m_l; m_f; theta_l; theta_f; d_l; d_f; r_f; xi_h; xi_l; k_h; k_l];
p = [m; g; l_0; m_l; m_f; d_l; d_f; r_f; xi_l; k_l];

%% Automatic Function Generation
% Generate MATLAB functions and export them to files
matlabFunction(f, 'File', 'dynamics', 'Vars', {x, u, p});
matlabFunction(M, 'File', 'massMatrix', 'Vars', {x, p});
matlabFunction(-C-G+Fspring_l+B*u, 'File', 'forces', 'Vars', {x, u, p});
matlabFunction(E, 'File', 'energyLevel', 'Vars', {x, u, p});
matlabFunction(posMechWork, 'File', 'posMechWork', 'Vars', {x, u});
matlabFunction(A_lin, B_lin, D_lin, 'File', 'get_linearDynamics', 'Vars', {p});

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

%% Define a function that indicates the position of the model states in the class states.
idx_x = [5; 10];
idx_u = [2];  % Example input indices

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
    C = simplify(expand(d_dTdqdt_dq * dq - dT_dq));  % Coriolis vector
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
