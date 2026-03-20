function lambda = compute_lambda(X_m, u, p)
% Calculate the contact force between the ground and the stance foot
import BF_run.SingleStance.*

% Get system dimensions
[~, nq_m, ~, ~] = get_sizes();

% Extract state from vector X
q_m = X_m(1:nq_m);
dq_m = X_m(nq_m+1:end);

% Extend to floating base coordinate (add x, y translation)
q_f = vertcat(0, 0, q_m);
dq_f = vertcat(0, 0, dq_m);
X_f = [q_f;dq_f];
nq_f = numel(q_f);

% Floating base mass matirx and forces
M_f = massMatrix_f(X_f, p);
F_f= forces_f(X_f, u, p);
% Floating base constraints jacobian
% Stance foot constraint
J_Stance = J_StanceConstraint_f(X_f, p);
dJ_Stance = dJ_StanceConstraint_f(X_f, p);

% [M, -J'][q_dot ] = [     F]
% [J,   0][lambda]   [-dJ*dq]
% gamma = dJ* dq
gamma = dJ_Stance * dq_f;
A = [M_f    , -J_Stance';
    J_Stance, zeros(size(J_Stance,1))];
B = [F_f;-gamma];

sol = A \ B;
lambda = sol(nq_f+1:end);
end
