function [X_post, lambda] = jump(X, u, p)
import BF_run.SingleStance.*;

% Extract from minimal coordinates
nq_m      = numel(X)/2;
X_m       = X;
q_pre_m     = X_m(1:nq_m);        %q
dq_pre_m    = X_m(nq_m+1:end);    %dq

% Extend to floating base coordinates
q_post_f  = [0;0;q_pre_m];
dq_post_f = [0;0;dq_pre_m];
X_post    = [q_post_f; dq_post_f];
% Compute lambda
lambda = compute_lambda(X, u, p);

end
