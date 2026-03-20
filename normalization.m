m   = 37;    % Total mass [kg]
g   = 9.81;  % Gravity [m/s²]
l   = 0.8;   % Leg length [m]

v = 0.6;
scaler_t = sqrt(l/g);
v_normalized1 = v*scaler_t/l;
disp(v_normalized);
