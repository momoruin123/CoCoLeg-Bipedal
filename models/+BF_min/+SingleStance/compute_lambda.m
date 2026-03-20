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

%% function lambda = compute_lambda(X_m, u, p)
% import BF_min.SingleStance.*;
%
% %Calculate the contact force between the ground and the stance foot
%
% [nx, nq, nu, np] = get_sizes();
%
% % Extract state from vector X
% q_m = X_m(1:nq);
% dq_m = X_m(nq+1:end);
%
% % Determine if it is an impact segment
% [swingFoot_p, swingFoot_v] = get_swingFoot_state(X_m, u, p);
%
% % Extend to floating base coordinate
% q_f = vertcat(0, 0, q_m);
% dq_f = vertcat(0, 0, dq_m);
% nq_f = numel(q_f);
% X_f = [q_f;dq_f];
%
% M_f = massMatrix_f(X_f, p);
% F_f= forces_f(X_f, u, p);
%
% if (swingFoot_p(2) < 1e-8) && (swingFoot_v(2) < 0)
%     J = J_SwingConstraint_f(X_f, p);
%     dJ = dJ_SwingConstraint_f(X_f, p);
% else
%     J = J_StanceConstraint_f(X_f, p);
%     dJ = dJ_StanceConstraint_f(X_f, p);
% end
%
% % gamma = dJ* dq
% gamma = dJ * dq_f;
%
% A = [M_f    , -J';
%     J, zeros(size(J,1))];
% B = [F_f;
%     -gamma];
%
% sol = A \ B;
% lambda = sol(nq_f+1:end);
%
%     function [SwingFoot_p, SwingFoot_v] = get_swingFoot_state(X_m, u, p)
%
%         % Extend to minimal coordinate
%         nq_m      = numel(X_m)/2;
%         q_pre_m     = X_m(1:nq_m);        %q
%         dq_pre_m    = X_m(nq_m+1:end);    %dq
%
%         % Compute swing foot position and velocity
%         J_swing  = J_SwingConstraint_m(X_m, p);
%         SwingFoot_p = cons_SwingFoot_m(q_pre_m, p);
%         SwingFoot_v = J_swing*dq_pre_m;
%
%     end
%
% end


%% function lambda = compute_lambda(X_m, u, p)
% import BF_min.SingleStance.*;
%
% %Calculate the contact force between the ground and the stance foot
%
% [nx, nq, nu, np] = get_sizes();
%
% % Extract state from vector X
% q_m = X_m(1:nq);
% dq_m = X_m(nq+1:end);
%
% % Extend to floating base coordinate
% extend_state = [0; 0];
% q_f = vertcat(extend_state, q_m);
% dq_f = vertcat(extend_state, dq_m);
% nq_f = numel(q_f);
% X_f = [q_f;dq_f];
%
% M_f = massMatrix_f(X_f, p);
% F_f= forces_f(X_f, u, p);
%
% J_Stance = J_StanceConstraint_f(X_f, p);
% dJ_Stance = dJ_StanceConstraint_f(X_f, p);
%
% % gamma = dJ* dq
% gamma = dJ_Stance * dq_f;
%
% A = [M_f    , -J_Stance';
%     J_Stance, zeros(size(J_Stance,1))];
% B = [F_f;
%     -gamma];
%
% sol = A \ B;
% lambda = sol(nq_f+1:end);
% end