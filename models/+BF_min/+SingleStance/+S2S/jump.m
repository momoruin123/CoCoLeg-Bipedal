function [X_post, lambda] = jump(X, u, p)
import BF_min.SingleStance.*;

% Extract from minimal coordinate
nq_m      = numel(X)/2;
X_m       = X;
q_pre_m     = X_m(1:nq_m);        %q
dq_pre_m    = X_m(nq_m+1:end);    %dq

% Extend to floating base coordinate and switch stance foot
T = [1, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1;
    0, 0, 0, 0, 0, 1, 0;];

q_pre_f     = [0;0;q_pre_m];
dq_pre_f    = [0;0;dq_pre_m];
X_f         = [q_pre_f;dq_pre_f];
nq_f        = numel(q_pre_f);

% compute q_post
q_post_f = T * q_pre_f;

% compute dq_post
J_swing = J_SwingConstraint_f(X_f, p); % constraint kinematics
J = J_swing;
nc = size(J, 1);

M = massMatrix_f(X_f, p); % mass matrix
A = [M, -J'; J, zeros(nc)];
B = [M*dq_pre_f; zeros(nc, 1)];
sol = A \ B;

dq_post_f = T * sol(1:nq_f);
lambda = sol(nq_f+1:end);
min_index = nq_f-nq_m+1;
q_post_m  = q_post_f(min_index:end);
dq_post_m = dq_post_f(min_index:end);
X_post = [q_post_m; dq_post_m];

end

%% function [X_post, lambda] = jump(X, u, p)
% % When event happens"
% % 1. Compute position and velocity of swing foot with swing foot jacobian;
% % 2. Swith stance foot. old swing foot as new stance foot, old stance foot
% %    as new swing foot
% % 3. Use stance foot jacobian compute acclucate and impluse after impact.
% 
% import BF_min.SingleStance.*;
% 
% % Extend to minimal coordinate
% nq_m      = numel(X)/2;
% X_m       = X;
% q_pre_m     = X_m(1:nq_m);        %q
% dq_pre_m    = X_m(nq_m+1:end);    %dq
% 
% % Compute swing foot position and velocity
% J_swing  = J_SwingConstraint_m(X_m, p); 
% SwingFoot_v = J_swing*dq_pre_m;
% SwingFoot_p = cons_SwingFoot_m(q_pre_m, p);
% 
% % Extend to floating base coordinate and switch stance foot
% T = [1, 0, 0, 0, 0;
%     0, 0, 1, 0, 0;
%     0, 1, 0, 0, 0;
%     0, 0, 0, 0, 1;
%     0, 0, 0, 1, 0;];
% q_pre_f     = [SwingFoot_p(1);SwingFoot_p(2);T*q_pre_m];
% dq_pre_f    = [SwingFoot_v(1);SwingFoot_v(2);T*dq_pre_m];
% X_f         = [q_pre_f;dq_pre_f];
% nq_f        = numel(q_pre_f);
% 
% % compute q_post
% q_post_f = q_pre_f;
% 
% % compute dq_post
% J_stance = J_StanceConstraint_f(X_f, p); % constraint kinematics
% J = J_stance;
% nc = size(J, 1);
% 
% M = massMatrix_f(X_f, p); % mass matrix
% A = [M, -J'; J, zeros(nc)];
% B = [M*dq_pre_f; zeros(nc, 1)];
% sol = A \ B;
% 
% dq_post_f = sol(1:nq_f);
% lambda = sol(nq_f+1:end);
% 
% min_index = nq_f-nq_m+1;
% q_post_m  = q_post_f(min_index:end);
% dq_post_m = dq_post_f(min_index:end);
% X_post = [q_post_m; dq_post_m];
% 
% end

%% function X_post = jump(X, u, p)
%     %
%     import BF_min.SingleStance.*;
%     nq        = numel(X)/2;
%     q_pre     = X(1:nq);
%     dq_pre    = X(nq+1:end);
%     q_post    = q_pre;
%
%     % Extend to floating base coordinate
%     q_pre_f   = [0;0;q_pre];
%     dq_pre_f  = [0;0;dq_pre];
%     X_f       = [q_pre_f;dq_pre_f];
%     nq_f      = numel(q_pre_f);
%     q_post_f  = q_pre_f;
%
%
%     T = [1, 0, 0, 0, 0, 0, 0;
%          0, 1, 0, 0, 0, 0, 0;
%          0, 0, 1, 0, 0, 0, 0;
%          0, 0, 0, 0, 1, 0, 0;
%          0, 0, 0, 1, 0, 0, 0;
%          0, 0, 0, 0, 0, 0, 1;
%          0, 0, 0, 0, 0, 1, 0;];
%
%     % compute dq_post
%     J_swing  = J_SwingConstraint_f(X_f, p); % constraint kinematics
%     J_stance = J_StanceConstraint_f(X_f, p);
%     J = J_swing;
%     nc = size(J, 1);
%
%     M = massMatrix_f(X_f, p); % mass matrix
%     A = [M, -J'; J, zeros(nc)];
%     B = [M*dq_pre_f; zeros(nc, 1)];
%     sol = A \ B;
%
%     dq_post_f = sol(1:nq_f);
%     lambda = sol(nq_f+1:end);
%
%     % Exchange left and right leg
%     FootPos = FootPositions(X_f,p);a
%     SwingFoot = FootPos(:,2);
%     q_post_f(1) = SwingFoot(1);
%
%
%     q_post_f  = T*q_post_f;
%     dq_post_f = T*dq_post_f;
%
%     q_post  = q_post_f(3:end);
%     dq_post = dq_post_f(3:end);
%     X_post = [q_post; dq_post];
%     % v_SwingFoot = J_swing * dq_post;
%     % FootPos = FootPtPositions(X_post,p);
%     % v_StanceFoot = J_stance * dq_post
% end