function f_S = dynamics(X, u, p)
import BF_run.SingleStance.*;

% Initialization
nq = length(X)/2;
q = X(1:nq);
dq = X(nq+1:end);

% ddq
% cumpute ddq under minimal coordinate
M = massMatrix_m(X, p);
F = forces_m(X, u, p);
ddq = M \ F;

% dynamic
f_S = [dq;ddq];
end
