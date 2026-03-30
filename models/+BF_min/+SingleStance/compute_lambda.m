%% Calculate the contact force between the ground and the stance foot
function lambda = compute_lambda(X_m, u, p)
import BF_min.SingleStance.*

[~, nq_m, ~, ~] = get_sizes();

% Extract state from vector X
q_m = X_m(1:nq_m);
dq_m = X_m(nq_m+1:end);

temp = [0; 0];
q_f = vertcat(temp, q_m);
dq_f = vertcat(temp, dq_m);
nq_f = numel(q_f);
X_f = [q_f;dq_f];

M_f = massMatrix_f(X_f, p);
F_f= forces_f(X_f, u, p);

J_Stance = J_StanceConstraint_f(X_f, p);
dJ_Stance = dJ_StanceConstraint_f(X_f, p);

% gamma = dJ* dq
gamma = dJ_Stance * dq_f;

A = [M_f    , -J_Stance';
    J_Stance, zeros(size(J_Stance,1))];
B = [F_f;
    -gamma];

sol = A \ B;
lambda = sol(nq_f+1:end);
end
