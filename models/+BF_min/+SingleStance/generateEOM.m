% System: bipedal robot with only SingleStance phases
% phase: SingleStance
% Generates dynamics, constraints, and transition functions for optimization

cd(fileparts(mfilename('fullpath')))
%% Floating-Base coordinate description
% Generalized coordinates and velocities
x       = sym('x', 'real');           % x of 
y       = sym('y', 'real');           % y of 
phi     = sym('phi', 'real');         % Pitch of Torso
alpha_l = sym('alpha_l','real');      % Hip_leftLeg angle
alpha_r = sym('alpha_r','real');      % Hip_rightLeg angle
beta_l  = sym('beta_l','real');       % leftLeg_knee angle
beta_r  = sym('beta_r','real');       % rightLeg_knee angle

dx       = sym('dx', 'real');         % dx of 
dy       = sym('dy', 'real');         % dy of 
dphi     = sym('dphi', 'real');       % Pitch of Torso
dalpha_l = sym('dalpha_l','real');    % Hip_leftLeg angle
dalpha_r = sym('dalpha_r','real');    % Hip_rightLeg angle
dbeta_l  = sym('dbeta_l','real');     % leftLeg_knee angle
dbeta_r  = sym('dbeta_r','real');     % rightLeg_knee angle

% State vectors
q_f  = [x y phi alpha_l alpha_r beta_l beta_r]';        
dq_f = [dx dy dphi dalpha_l dalpha_r dbeta_l dbeta_r]';   

%% Minimal coordinate description
q1 = phi;
q2 = alpha_l;
q3 = alpha_r;
q4 = beta_l;
q5 = beta_r;

dq1 = dphi;
dq2 = dalpha_l;
dq3 = dalpha_r;
dq4 = dbeta_l;
dq5 = dbeta_r;

% State vectors
q_m  = [q1 q2 q3 q4 q5]';        % Generalized coordinates
dq_m = [dq1 dq2 dq3 dq4 dq5]';   % Generalized velocities

%% Control inputs
u_hl = sym('u_hl','real');  % Hip_leftLeg torque
u_hr = sym('u_hr','real');  % Hip_rightLeg torque
u_kl = sym('u_kl','real');   % LeftLeg_knee torque
u_kr = sym('u_kr','real');   % RightLeg_knee torque

%% State vectors
u = [u_hl; u_hr; u_kl; u_kr]; % input vector

X_f = [q_f; dq_f];
X_m = [q_m; dq_m];   

%% Model Parameters
% Define symbolic physical parameters
g       = sym('g','real');          % Gravity
m       = sym('m','real');          % Total mass
l       = sym('l','real');          % leg length

m_2     = sym('m_2','real');        % Thigh mass
m_3     = sym('m_3','real');        % Shank mass
theta_1 = sym('theta_1','real');    % Torso inertia
theta_2 = sym('theta_2','real');    % Thigh inertia
theta_3 = sym('theta_3','real');    % Shank inertia
l_1     = sym('l_1','real');        % Torso length
l_2     = sym('l_2','real');        % Thigh length
d_1     = sym('d_1','real');        % distance from torso CoM to hip
d_2     = sym('d_2','real');        % distance from hip to thigh CoM
d_3     = sym('d_3','real');        % distance from knee to shank CoM

% Spring parameters
xi_h    = sym('xi_h','real');       % Hip damping ratio
xi_k    = sym('xi_k','real');       % Knee damping ratio
k_h     = sym('k_h','real');        % Hip stiffness
k_k     = sym('k_k','real');        % Knee stiffness
alpha0  = sym('alpha0','real');     % Hip spring reference position
beta0   = sym('beta0','real');      % Kenn spring reference position

% Indirect variables
m_1     = m - 2*m_2 - 2*m_3;        % Torso mass
l_3     = l - l_2;                  % Shank length

% Gravity vector
g_vec   = [0; -g];

%% DYNAMICS (obtained via the Euler-Lagrange equation)
% CoG-orientations (from kinematics):
CoG_Torso_ang       = q1;
CoG_StanceThigh_ang = q1 + q2;
CoG_StanceShank_ang = q1 + q2 + q4;
CoG_SwingThigh_ang  = q1 + q3;
CoG_SwingShank_ang  = q1 + q3 + q5;

% Links-positions (from kinematics):
StanceFoot  = [x; y];
StanceKnee  = StanceFoot + [-l_3*sin(CoG_StanceShank_ang);
                            +l_3*cos(CoG_StanceShank_ang)];
Hip         = StanceKnee + [-l_2*sin(CoG_StanceThigh_ang);
                            +l_2*cos(CoG_StanceThigh_ang)];
Head        = Hip + [-l_1*sin(CoG_Torso_ang);
                     +l_1*cos(CoG_Torso_ang)];

SwingKnee   = Hip + [+l_2*sin(CoG_SwingThigh_ang);
                     -l_2*cos(CoG_SwingThigh_ang)];
SwingFoot   = SwingKnee + [+l_3*sin(CoG_SwingShank_ang);
                           -l_3*cos(CoG_SwingShank_ang)];

% CoG-positions (from kinematics):
CoG_StanceShank = StanceFoot + [-(l_3-d_3)*sin(CoG_StanceShank_ang);
                                +(l_3-d_3)*cos(CoG_StanceShank_ang)];
CoG_StanceThigh = StanceKnee + [-(l_2-d_2)*sin(CoG_StanceThigh_ang);
                                +(l_2-d_2)*cos(CoG_StanceThigh_ang)];
CoG_Torso       = Hip + [-d_1*sin(CoG_Torso_ang);
                         +d_1*cos(CoG_Torso_ang)];

CoG_SwingThigh  = Hip + [+d_2*sin(CoG_SwingThigh_ang);
                         -d_2*cos(CoG_SwingThigh_ang)];
CoG_SwingShank  = SwingKnee + [+d_3*sin(CoG_SwingShank_ang);
                               -d_3*cos(CoG_SwingShank_ang)];

% CoG-velocities (computed via jacobians):
d_CoG_Torso           = jacobian(CoG_Torso,q_f)*dq_f;
d_CoG_StanceThigh     = jacobian(CoG_StanceThigh,q_f)*dq_f;
d_CoG_StanceShank     = jacobian(CoG_StanceShank,q_f)*dq_f;
d_CoG_SwingThigh      = jacobian(CoG_SwingThigh,q_f)*dq_f;
d_CoG_SwingShank      = jacobian(CoG_SwingShank,q_f)*dq_f;
d_CoG_SwingFoot       = jacobian(SwingFoot,q_f)*dq_f;

d_CoG_Torso_ang       = jacobian(CoG_Torso_ang,q_f)*dq_f;
d_CoG_StanceThigh_ang = jacobian(CoG_StanceThigh_ang,q_f)*dq_f;
d_CoG_StanceShank_ang = jacobian(CoG_StanceShank_ang,q_f)*dq_f;
d_CoG_SwingThigh_ang  = jacobian(CoG_SwingThigh_ang,q_f)*dq_f;
d_CoG_SwingShank_ang  = jacobian(CoG_SwingShank_ang,q_f)*dq_f;

% Distance from Hip to CoG_Shank
Hip2StanceShank    = Hip-CoG_StanceShank;
Hip2StanceShank_sq = simplify(Hip2StanceShank'*Hip2StanceShank);
Hip2SwingShank     = Hip-CoG_SwingShank;
Hip2SwingShank_sq  = simplify(Hip2SwingShank'*Hip2SwingShank);    

% Damping coefficients from damping ratios
% b_hL = 2*xi_h*sqrt(k_h*(theta_2+m_2*d_2^2+theta_3+m_3*Hip2StanceShank_sq));
% b_hR = 2*xi_h*sqrt(k_h*(theta_2+m_2*d_2^2+theta_3+m_3*Hip2SwingShank_sq));
% b_kL = 2*xi_k*sqrt(k_k*(theta_3+m_3*d_3^2));
% b_kR = 2*xi_k*sqrt(k_k*(theta_3+m_3*d_3^2));

b_hL = 16.5;
b_hR = 16.5;
b_kL = 5.48;
b_kR = 5.48;

% b_hL = 0;
% b_hR = 0;
% b_kL = 0;
% b_kR = 0;

%% Summary of positons
% CoGs (incl orientation of the segments):
CoGs = [CoG_Torso,     CoG_StanceThigh,     CoG_StanceShank,     CoG_SwingThigh,     CoG_SwingShank;
        CoG_Torso_ang, CoG_StanceThigh_ang, CoG_StanceShank_ang, CoG_SwingThigh_ang, CoG_SwingShank_ang];

% Links (or rather: joint positions)
links = [StanceFoot, StanceKnee, Hip, Head, SwingKnee, SwingFoot];

% Position of the foot points:     
footPts = [StanceFoot, SwingFoot];

%% Energy Formulation
% Potential Energy (due to gravity):
V_grav = CoG_Torso(2)*m_1*g + ...
    CoG_StanceThigh(2)*m_2*g + ...
    CoG_SwingThigh(2)*m_2*g + ...
    CoG_StanceShank(2)*m_3*g + ...
    CoG_SwingShank(2)*m_3*g;
V_grav = simplify(V_grav);
V_spring  = 0.5*k_h*(alpha0-alpha_l).^2 + ...
            0.5*k_h*(alpha0-alpha_r).^2 + ...
            0.5*k_k*(beta0-beta_l).^2 + ...
            0.5*k_k*(beta0-beta_r).^2;
V = V_grav + V_spring;

% Kinetic Energy:
T = 0.5 * ( m_1 * (d_CoG_Torso'*d_CoG_Torso) + ...
            m_2 * (d_CoG_StanceThigh'*d_CoG_StanceThigh) + ...
            m_2 * (d_CoG_SwingThigh'*d_CoG_SwingThigh) + ...
            m_3 * (d_CoG_StanceShank'*d_CoG_StanceShank) + ...
            m_3 * (d_CoG_SwingShank'*d_CoG_SwingShank) + ...
            theta_1 * d_CoG_Torso_ang^2 + ...
            theta_2 * d_CoG_StanceThigh_ang^2 + ...
            theta_2 * d_CoG_SwingThigh_ang^2 + ...
            theta_3 * d_CoG_StanceShank_ang^2 + ...
            theta_3 * d_CoG_SwingShank_ang^2);
T = simplify(T);

%% Total Energy
E = simplify(T + V);

%% Euler-Lagrange Equations
% M * ddq + C * dq + G = 0
% Floating base 
[M_f, C_f, G_f] = eulerLagrange(T, V_grav, q_f, dq_f); 

% Mininal
T = simplify(subs(T, {x, y, dx, dy}, {0, 0, 0, 0}));
V_grav = simplify(subs(V_grav, {x, y, dx, dy}, {0, 0, 0, 0}));
[M_m, C_m, G_m] = eulerLagrange(T, V_grav, q_m, dq_m); 

% Spring and damping forces
Fspring_leftHip   = k_h*(alpha0-alpha_l) + b_hL*(0-dalpha_l); 
Fspring_rightHip  = k_h*(alpha0-alpha_r) + b_hR*(0-dalpha_r); 
Fspring_leftKnee  = k_k*(beta0-beta_l)  + b_kL*(0-dbeta_l); 
Fspring_rightKnee = k_k*(beta0-beta_r)  + b_kR*(0-dbeta_r); 

% Input mapping matrices
% Floating base input matrix
B_f = [[0,0,0,0];...      % x
       [0,0,0,0];...      % y
       [0,0,0,0];...      % phi
       [1,0,0,0];...      % alpha_l
       [0,1,0,0];...      % alpha_r
       [0,0,1,0];...      % beta_l
       [0,0,0,1]];        % beta_r
% Minimal coordinate input matrix
B_m = [[0,0,0,0];...
       [1,0,0,0];...
       [0,1,0,0];...
       [0,0,1,0];...
       [0,0,0,1]];

% Passive spring torques
tau_Spring_f = B_f*[Fspring_leftHip; Fspring_rightHip; Fspring_leftKnee; Fspring_rightKnee];
tau_Spring_m = B_m*[Fspring_leftHip; Fspring_rightHip; Fspring_leftKnee; Fspring_rightKnee];

%% Continuous Dynamics
% Compute acceleration: ddq = M^(-1) * (tau - C*dq - G)
% acc = M \ (-C-G+tau_J+B*u);
forces_f = simplify(-C_f-G_f + tau_Spring_f + B_f*u);
forces_m = simplify(-C_m-G_m + tau_Spring_m + B_m*u);

%% CONSTRAINT DYNAMICS 
% Contact points:
pos_StanceFoot_f = StanceFoot;
pos_SwingFoot_f  = SwingFoot;
pos_StanceFoot_m = subs(pos_StanceFoot_f, {x, y, dx, dy}, {0, 0, 0, 0});
pos_SwingFoot_m  = subs(pos_SwingFoot_f, {x, y, dx, dy}, {0, 0, 0, 0});

% Constraint Jacobians:
% Floating base
[J_StanceConstraint_f, dJ_StanceConstraint_f] = computeContactMatrix(pos_StanceFoot_f, q_f, dq_f);
[J_SwingConstraint_f , dJ_SwingConstraint_f] = computeContactMatrix(pos_SwingFoot_f, q_f, dq_f);
% Minimal coordinate
[J_StanceConstraint_m, dJ_StanceConstraint_m] = computeContactMatrix(pos_StanceFoot_m, q_m, dq_m);
[J_SwingConstraint_m , dJ_SwingConstraint_m] = computeContactMatrix(pos_SwingFoot_m, q_m, dq_m);

%% Inverse Kinematics
Hip_m = subs(Hip, {x, y, dx, dy}, {0, 0, 0, 0});
[J_Hip_m, ~] = computeContactMatrix(Hip_m, q_m, dq_m);

%% Event and Jump Functions for Phase Transitions
event_td = SwingFoot(2);
event_td = subs(event_td, {x, y, dx, dy}, {0, 0, 0, 0});

%% Compute positive mechanical work
% Smooth approximation of mechanical power
posMechWork = smoothRampFcn(u'*dq_m(2:end), 0.3);

%% parameters vector
p = [g; m; l; m_2; m_3; theta_1; theta_2; theta_3; l_1; l_2; d_1; d_2; d_3; alpha0; beta0; xi_h; xi_k; k_h; k_k];  % k_h and k_k are the only params to be optimized on here

%% Automatic Function Generation
% Generate MATLAB functions for the dynamics and their derivatives
% Generate MATLAB functions and export them to files
matlabFunction(M_f, 'File', 'massMatrix_f', 'Vars', {X_f, p});
matlabFunction(forces_f, 'File', 'forces_f', 'Vars', {X_f, u, p});
matlabFunction(M_m, 'File', 'massMatrix_m', 'Vars', {X_m, p});
matlabFunction(forces_m, 'File', 'forces_m', 'Vars', {X_m, u, p});

matlabFunction(E, 'File', 'energyLevel', 'Vars', {X_f, u, p});
matlabFunction(posMechWork, 'File', 'posMechWork', 'Vars', {X_m, u});

% for constraints and collisions:
matlabFunction(pos_StanceFoot_f,'file','pos_StanceFoot_f','vars',{X_f, p});
matlabFunction(J_StanceConstraint_f,'file','J_StanceConstraint_f','vars',{X_f, p});
matlabFunction(dJ_StanceConstraint_f,'file','dJ_StanceConstraint_f','vars',{X_f p});
matlabFunction(pos_StanceFoot_m,'file','pos_StanceFoot_m','vars',{X_m, p});
matlabFunction(J_StanceConstraint_m,'file','J_StanceConstraint_m','vars',{X_m, p});
matlabFunction(dJ_StanceConstraint_m,'file','dJ_StanceConstraint_m','vars',{X_m p});

matlabFunction(pos_SwingFoot_f,'file','pos_SwingFoot_f','vars',{X_f, p});
matlabFunction(J_SwingConstraint_f,'file','J_SwingConstraint_f','vars',{X_f, p});
matlabFunction(dJ_SwingConstraint_f,'file','dJ_SwingConstraint_f','vars',{X_f p});
matlabFunction(pos_SwingFoot_m,'file','pos_SwingFoot_m','vars',{X_m, p});
matlabFunction(J_SwingConstraint_m,'file','J_SwingConstraint_m','vars',{X_m, p});
matlabFunction(dJ_SwingConstraint_m,'file','dJ_SwingConstraint_m','vars',{X_m p});

% for inverse kinematics
matlabFunction(J_Hip_m,'file','J_Hip_m','vars',{X_m, p});

% for graphics:
matlabFunction(CoGs,'file','CoGPositions_f','vars',{X_f, p});
matlabFunction(links,'file','LinkPositions_f','vars',{X_f, p});
matlabFunction(footPts,'file','FootPositions_f','vars',{X_f, p});

% spring forces
matlabFunction(tau_Spring_f,'file','tau_Spring_f','vars',{X_f, p});
matlabFunction(tau_Spring_m,'file','tau_Spring_m','vars',{X_m, p});

%% Create packages for transitions
% Define package folder names (in current folder)
packageNames = {'S2S'};

% Create package folders in current directory if they don't exist
for i = 1:length(packageNames)
    packageName = packageNames{i};
    packageFolder = ['+', packageName];

    if ~exist(packageFolder, 'dir')
        mkdir(packageFolder);
    end
end

% Generate the MATLAB functions in the appropriate package folders
matlabFunction(event_td, 'File', fullfile(['+', packageNames{1}], 'event'), 'Vars', {X_m, u, p});

%% Define sizes in a function
nX_m = numel(X_m);
nq_m = numel(q_m);
nu   = numel(u);
np   = numel(p);

% Create a function file that returns these sizes
fid = fopen('get_sizes.m', 'w');
fprintf(fid, 'function [nX_m, nq_m, nu, np] = get_sizes()\n');
fprintf(fid, 'nX_m = %d;\n', nX_m);
fprintf(fid, 'nq_m = %d;\n', nq_m);
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

%% Define a function that indicates the position of the model states in the class states.
idx_x = [5; 10];
idx_u = [4];  % Example input indices

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
    % using the Euler-Lagra nge formalism.

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

function [J, dJ] = computeContactMatrix(con, q, dq)
    % Compute the contact projection matrix W and its time derivative W_dot
    % for a given constraint 'con'.
    
    J = jacobian(con, q)';  % Contact matrix
    dJ = J;

    % Compute time derivative of each column of W
    for j = 1:size(J, 2)
        dJ(:, j) = jacobian(J(:, j), q) * dq;  % Time derivative of contact matrix
    end
    J = J';
    dJ = dJ';
end
