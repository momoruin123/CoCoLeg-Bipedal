function event_lo = event(X, u, p)
import BF_run.SingleStance.*
% lift off event
% Compute ground reaction force and treat it as event 
lambda = compute_lambda(X, u, p);
event_lo = lambda(2);
end
