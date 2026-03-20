function f_S = dynamics(X_f, u, p)
% Now the system is in floating base coordinates
import BF_run.SingleStance.*;

% Initialization
nq = length(X_f)/2;
q = X_f(1:nq);
dq = X_f(nq+1:end);

% ddq
M = massMatrix_f(X_f, p);
F = forces_f(X_f, u, p);
ddq = M \ F;

% dynamic
f_S = [dq;ddq];
end
