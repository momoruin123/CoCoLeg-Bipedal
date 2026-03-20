function [eventFcn, jumpMap, nX_post] = createEventJump()
import casadi.*
import BF_run.SingleStance.*;

% Get system dimensions
[nX, nq, nu, np] = get_sizes();
states     = MX.sym('states', nX);
inputs     = MX.sym('inputs', nu);
parameters = MX.sym('parameters', np);

% Original minimal coordinate dimensions
q_pre_m    = states(1:nq);
dq_pre_m   = states(nq+1:end);
X_m        = [q_pre_m; dq_pre_m];

% Extend to floating base coordinate (add x, y translation)
T = [1, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1;
    0, 0, 0, 0, 0, 1, 0;];

q_pre_f   = vertcat(0, 0, q_pre_m);   % [x; y; q_joint]
dq_pre_f  = vertcat(0, 0, dq_pre_m);  % [dx; dy; dq_joint]
X_f       = vertcat(q_pre_f, dq_pre_f);
nq_f      = numel(q_pre_f);

% compute q_post
q_post_f = T * q_pre_f;

% compute dq_post
% Floating base mass matirx and constraint jacobian matrix
M  = massMatrix_f(X_f, parameters);
J_swing = J_SwingConstraint_f(X_f, parameters);
J = J_swing;
nc = size(J, 1);

% Contact equation
A = vertcat( horzcat(M, -J'), ...
    horzcat(J, MX.zeros(nc, nc)) );
B = vertcat( M*dq_pre_f, ...
    MX.zeros(nc, 1) );
plus = solve(A, B);

% dq_post_f
dq_post_f = T * plus(1:nq_f);
% dq_post_f = plus(1:nq_f);

% Remove the floating base and return to the state vector
% used for optimization.
q_post  = q_post_f(3:end);
dq_post = dq_post_f(3:end);
g = [q_post; dq_post];

nX_post = numel(g);
e = S2S.event(states, inputs, parameters);  % Touchdown event

% Create CasADi functions
eventFcn = Function('eventFcn', {states, inputs, parameters}, {e});
jumpMap  = Function('jumpMap', {states, inputs, parameters}, {g});
end

%% function [eventFcn, jumpMap, nX_post] = createEventJump(p)
% import casadi.*
% import BF_min.SingleStance.*;
% 
% % Get system dimensions
% [nX, nq, nu, np] = get_sizes();
% states     = MX.sym('states', nX);
% inputs     = MX.sym('inputs', nu);
% parameters = MX.sym('parameters', np);
% 
% % Original minimal coordinate dimensions
% q_pre_m    = states(1:nq);
% dq_pre_m   = states(nq+1:end);
% X_m        = [q_pre_m; dq_pre_m];
% 
% % Compute swing foot position and velocity
% J_swing  = J_SwingConstraint_m(X_m, parameters); 
% SwingFoot_v = J_swing*dq_pre_m;
% SwingFoot_p = cons_SwingFoot_m(q_pre_m, parameters);
% 
% % Extend to floating base coordinate (add x, y translation)
% T = [1, 0, 0, 0, 0;
%     0, 0, 1, 0, 0;
%     0, 1, 0, 0, 0;
%     0, 0, 0, 0, 1;
%     0, 0, 0, 1, 0;];
% 
% q_pre_f   = vertcat(SwingFoot_p(1), SwingFoot_p(2), T*q_pre_m);   % [x; y; q_joint]
% dq_pre_f  = vertcat(SwingFoot_v(1), SwingFoot_v(2), T*dq_pre_m);  % [dx; dy; dq_joint]
% X_f       = vertcat(q_pre_f, dq_pre_f);
% nq_f      = numel(q_pre_f);
% 
% % compute q_post
% q_post_f = q_pre_f;
% 
% % compute dq_post
% % Floating base mass matirx and constraint jacobian matrix
% M  = massMatrix_f(X_f, parameters);
% J_stance = J_StanceConstraint_f(X_f, parameters);
% J = J_stance;
% nc = size(J, 1);
% % Contact equation
% A = vertcat( horzcat(M, -J'), ...
%     horzcat(J, MX.zeros(nc, nc)) );
% B = vertcat( M*dq_pre_f, ...
%     MX.zeros(nc, 1) );
% plus = solve(A, B);
% % dq_post_f
% dq_post_f = plus(1:nq_f);
% 
% % Remove the floating base and return to the state vector
% % used for optimization.
% q_post  = q_post_f(3:end);
% dq_post = dq_post_f(3:end);
% g = [q_post; dq_post];
% 
% nX_post = numel(g);
% e = S2S.event(states, inputs, parameters);  % Touchdown event
% 
% % Create CasADi functions
% eventFcn = Function('eventFcn', {states, inputs, parameters}, {e});
% jumpMap  = Function('jumpMap', {states, inputs, parameters}, {g});
% end

%%
% function [eventFcn, jumpMap, nX_post] = createEventJump()
% import casadi.*
% import BF_min.SingleStance.*;
% 
% % Get system dimensions
% [nX, nq, nu, np] = get_sizes();
% states     = MX.sym('states', nX);
% inputs     = MX.sym('inputs', nu);
% parameters = MX.sym('parameters', np);
% 
% % Original minimal coordinate dimensions
% q_pre     = states(1:nq);
% dq_pre    = states(nq+1:end);
% 
% % Extend to floating base coordinate (add x, y translation)
% q_pre_f   = vertcat(0, 0, q_pre);   % [x; y; q_joint]
% dq_pre_f  = vertcat(0, 0, dq_pre);  % [dx; dy; dq_joint]
% X_f       = vertcat(q_pre_f, dq_pre_f);
% nq_f      = numel(q_pre_f);
% 
% % Floating base mass matirx and constraint jacobian matrix
% M  = massMatrix_f(X_f, parameters);
% J_stance = J_StanceConstraint_f(X_f, parameters);
% J_swing  = J_SwingConstraint_f(X_f, parameters);
% % J = horzcat(J_stance, J_swing(:,2)); % Only retain y direction for swing foot
% J = J_swing;
% nc = size(J, 1);
% 
% % Contact equation
% A = vertcat( horzcat(M, -J'), ...
%     horzcat(J, MX.zeros(nc, nc)) );
% B = vertcat( M*dq_pre_f, ...
%     MX.zeros(nc, 1) );
% sol = A \ B;
% dq_post_f = sol(1:nq_f);
% 
% % Switch legs
% FootPos = FootPositions(X_f, parameters);
% SwingFoot = FootPos(:, 2);
% q_post_f = q_pre_f;
% q_post_f(1) = SwingFoot(1); % Update base x position
% 
% T = MX([1, 0, 0, 0, 0;
%     0, 0, 1, 0, 0;
%     0, 1, 0, 0, 0;
%     0, 0, 0, 0, 1;
%     0, 0, 0, 1, 0]);
% 
% q_post_f  = T * q_post_f;
% dq_post_f = T * dq_post_f;
% 
% % Remove the floating base and return to the state vector
% % used for optimization.
% q_post  = q_post_f(3:end);
% dq_post = dq_post_f(3:end);
% g = vertcat(q_post, dq_post);
% 
% nX_post = numel(g);
% e = S2S.event(states, inputs, parameters);  % Touchdown event
% 
% % Create CasADi functions
% eventFcn = Function('eventFcn', {states, inputs, parameters}, {e});
% jumpMap  = Function('jumpMap', {states, inputs, parameters}, {g});
% end

%%
% function [eventFcn, jumpMap, nx_post] = createEventJump()
%     import casadi.*
%     % Import all public members from mypackage
%     import BF_min.left_Stance.*;
%
%     % find a way to get nx, nq, nu and np.
%     [nx, nq, nu, np] = get_sizes();
%
%     states     = MX.sym('states', nx);
%     inputs     = MX.sym('inputs', nu);
%     parameters = MX.sym('parameters', np);
%
%     e = L2R.event(states, inputs, parameters);
%     g = L2R.jump(states, inputs, parameters);
%
%     nx_post = numel(g);
%
%     eventFcn = Function('event', {states, inputs, parameters}, {e});
%     jumpMap = Function('jumpMap', {states, inputs, parameters}, {g});
% end
