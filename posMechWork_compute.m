configfile_BF_min;
import BF_min.SingleStance.*
[~, nq_m, nu, ~] = get_sizes();
dt = 3e-2;
params = getModelParameters(config, 0, 0);
paramName = {'m', 'g'};
p_idx = getParameterIndex(paramName, config);
m = params(p_idx(1));
g = params(p_idx(2));

for i = 1:15
    trajPlot(i) = interpolate_Z2traj(config, Z_array(i,:)', [], 0, dt);

    % Extract 
    x = trajPlot(i).x;
    dq_act = x(:, 7:10);
    u = trajPlot(i).u;
    t = trajPlot(i).t;

    % Positive mechanic work
    dt = t(2)-t(1);
    P = u .* dq_act;
    W_pos = sum(max(P,0), 'all') * dt;

    % Hip position
    X_init_m = x(1,:);
    X_init_f = [0, 0, x(1,1:5), 0, 0, x(1,6:end)]';
    X_end_f  = [0, 0, x(end,1:5), 0, 0, x(end,6:end)]';
    links_init = LinkPositions(X_init_f, params);
    links_end  = LinkPositions(X_end_f, params);
    hip_init = links_init(:,3);
    hip_end  = links_end(:,3);

    % Distance
    dx = hip_end(1)-hip_init(1);

    % CoT
    CoT(i) = W_pos / (m * g * dx); 
end

disp(CoT)