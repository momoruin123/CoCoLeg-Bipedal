function f = smoothRampFcn(x, kappa) 
% SMOOTHRAMPCON Smooth approximation of ramp function max(0,x)
%
%   Inputs:
%     x     - input values
%     kappa - smoothing parameter (smaller = sharper transition)
%
%   Output:
%     f     - smoothed ramp function values
%
%   author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    x   = x(:);
    % idx = x/kappa>700;
    conv = (0.5*(-x)+sqrt((0.5*(-x)).^2+kappa^2));
    %f   = kappa*log(1+exp(x/kappa));
    f = 0.5*x+sqrt((0.5*x).^2+kappa^2) + 0.1*conv;
   
end