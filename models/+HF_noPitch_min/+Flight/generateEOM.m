% System: hopping robot with fixed main body pitch
% phase: Flight
% Generates dynamics, constraints, and transition functions for optimization
%
% Author: Maximilian Raff, Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025
cd(fileparts(mfilename('fullpath')))
%% Floating Base description
nx = 8;  % Number of states
nq = 4;  % Number of generalized coordinates

% Generalized coordinates
x      = sym('x','real');      % Horizontal position
y      = sym('y','real');      % Vertical position  
alpha  = sym('alpha','real');  % Hip angle
l      = sym('l','real');      % Leg length
% Generalized velocities
dx     = sym('dx','real');     % Horizontal velocity
dy     = sym('dy','real');     % Vertical velocity
dalpha = sym('dalpha','real'); % Hip angular velocity
dl     = sym('dl','real');     % Leg extension velocity
% time
t      = sym('t', 'real');

% Control inputs
syms u_h u_l real  % Hip torque and leg force

% State vectors
q  = [x; y; alpha; l];    % Generalized coordinates
dq = [dx; dy; dalpha; dl]; % Generalized velocities

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

%% Contact Constraints
% Define initial contact point parameters
syms x_c alpha_c l_c real  % Initial contact position and configuration

% Stance phase constraints (rolling contact without slipping)
con_S = [x + r_f*(phi_0 + alpha) + l*sin(phi_0 + alpha) - x_c - (alpha_c + phi_0)*r_f - l_c*sin(alpha_c + phi_0);  % Horizontal rolling constraint
         y - l*cos(phi_0 + alpha) - r_f];  % Vertical ground contact constraint

% Compute contact projection matrix and its time derivative
[W_S, W_S_dot] = computeContactMatrix(con_S, q, dq);  % Rolling contact kinematics

%% Kinematics
% Compute positions and velocities of all body segments
Angle_l   = alpha + phi_0;  % Leg angle
Angle_f   = alpha + phi_0;  % Foot angle

% Center of gravity positions
CoG_t     = [x; y];  % Torso CoG
CoG_l     = CoG_t + d_l * [sin(Angle_l);-cos(Angle_l)];  % Leg CoG
CoG_f     = CoG_t + (l-d_f) * [sin(Angle_f);-cos(Angle_f)];  % Foot CoG

% Velocities via Jacobians
d_Angle_l = jacobian(Angle_l, q) * dq;
d_Angle_f = jacobian(Angle_f, q) * dq;
d_CoG_t   = jacobian(CoG_t, q) * dq;  % Velocity of the torso
d_CoG_l   = jacobian(CoG_l, q) * dq;  % Velocity of the leg
d_CoG_f   = jacobian(CoG_f, q) * dq;  % Velocity of the foot

%% Energy Formulation
% Potential Energy (gravity + springs)
V_grav = CoG_t(2)*m_t*g + CoG_l(2)*m_l*g + CoG_f(2)*m_f*g;
V_spring  = 0.5*k_l*(l_0-l).^2 + 0.5*k_h*(0-alpha).^2;
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
B = [[0;0;1;0],...
     [0;0;0;1]];

tau_J = B*[Fspring_h; Fspring_l];

%% Continuous Dynamics
% Compute acceleration: ddq = M^(-1) * (tau - C*dq - G)
acc = M \(-C-G+tau_J+B*[u_h;u_l]);
% State derivative
f_F = [dq;acc];

%% Define event and mapping touch-dowm
% F2S: touchdown: foot touches the ground
event_td = con_S(2); % transition into stance
% Compute discrete map for collision dynamics
plus = [M -W_S; W_S' zeros(2)] \ [M*dq;0;0];
jump_td  = [alpha; l; plus(3:4)]; 

% F2F: apex event heighest point in flight
event_apex = dy;    
jump_apex  = [q; dq];

%% Compute positive mechanical work
% Smooth approximation of mechanical power
posMechWork = smoothRampFcn(dl*u_l + dalpha*u_h, 0.3);

%% Energy level
E = T + V;
%% parameters vector
% parameters to be optimized should always placed AT THE END OF THE VECTOR
p = [m; g; l_0; m_l; m_f; theta_l; theta_f; d_l; d_f; r_f; phi_0; xi_h; xi_l; k_h; k_l];

%% Automatic Function Generation
% Generate MATLAB functions for the dynamics and their derivatives
x = [q; dq];   % State vector
u = [u_h;u_l]; % input vector

% Generate MATLAB functions and export them to files
matlabFunction(f_F, 'File', 'dynamics', 'Vars', {x, u, p});
matlabFunction(M, 'File', 'massMatrix', 'Vars', {x, p});
matlabFunction(-C-G+tau_J+B*u, 'File', 'forces', 'Vars', {x, u, p});
matlabFunction(E, 'File', 'energyLevel', 'Vars', {x, u, p});
matlabFunction(posMechWork, 'File', 'posMechWork', 'Vars', {x, u});
matlabFunction(W_S, 'File', 'contactJacobian', 'Vars', {x, p});
matlabFunction(W_S_dot, 'File', 'contactJacobianDerivative', 'Vars', {x, p});

%% Create packages for transitions
% Define package folder names (in current folder)
packageNames = {'F2S', 'F2F'};

% Create package folders in current directory if they don't exist
for i = 1:length(packageNames)
    packageName = packageNames{i};
    packageFolder = ['+', packageName];

    if ~exist(packageFolder, 'dir')
        mkdir(packageFolder);
    end
end

% Generate the MATLAB functions in the appropriate package folders
matlabFunction(event_td, 'File', fullfile(['+', packageNames{1}], 'event'), 'Vars', {x, u, p});
matlabFunction(jump_td, 'File', fullfile(['+', packageNames{1}], 'jump'), 'Vars', {x, u, p});
matlabFunction(event_apex, 'File', fullfile(['+', packageNames{2}], 'event'), 'Vars', {x, u, p});
matlabFunction(jump_apex, 'File', fullfile(['+', packageNames{2}], 'jump'), 'Vars', {x, u, p});

%% Define sizes in a function
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

%% Define the function that indicates the position of the optimal states and inputs in the full states
idx_x = [1;2;4;5;6;7;9;10];
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
