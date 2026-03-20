function [X_post_m, lambda] = jump(X_pre_f, u, p)
% Jump map from floating base coordinates to minimal coordinates
% X_pre_f: 7*1;
% X_post_m: 5*1;
import BF_run.Flight.*;

nX_f = numel(X_pre_f);
nq_f = nX_f/2;

% Extract q and dq
q_pre_f  = X_pre_f(1:nq_f);
dq_pre_f = X_pre_f(nq_f+1:end);

% Extend to floating base coordinate and switch stance foot
T = [1, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1;
    0, 0, 0, 0, 0, 1, 0;];

% compute q_post
q_post_f = T * q_pre_f;

% compute dq_post
J_swing = J_SwingConstraint_f(X_pre_f, p); % constraint kinematics
J = J_swing;
nc = size(J, 1);

M = massMatrix_f(X_pre_f, p); % mass matrix
A = [M, -J'; J, zeros(nc)];
B = [M*dq_pre_f; zeros(nc, 1)];
sol = A \ B;

dq_post_f = T * sol(1:nq_f);
lambda = sol(nq_f+1:end);
q_post_m  = q_post_f(3:end);
dq_post_m = dq_post_f(3:end);
X_post_m = [q_post_m; dq_post_m];

end
