function x = Bipedal_IK(config, step_length, hip_height, phi, v_target)
% Bipedal inverse kinematic function 
%
%   It will rebuild bipedal robot states in minimal coordinates.
%
%   Inputs:
%     config      - system configuration structure
%     step_length - step length
%     hip_height  - hip height
%     phi         - pitch angle of main body
%     v_target    - target average velocity of hip 
%
%   Output:
%     x           - The initial state obtained from the corresponding 
%                   parameters. x = [q;dq]


import BF_min.SingleStance.*
[p, paraNames] = getFullParameters(config);

% Find idx of leg length 'l' and thigh length 'l1'
idx_l  = strcmp(paraNames, 'l'); 
idx_l_2 = strcmp(paraNames, 'l_2');

% Extract values
leg_length = p(idx_l);
thigh_length = p(idx_l_2);
shank_length = leg_length - thigh_length;

%% Standce Triangle
% Compute the angle of standce triangle, the vertex angle is 'B', the base 
% angle is 'A', the legs of isosceles triangle is 'a', the base is 'b'.

% Define the heigth and base of stance triangle
st_h = hip_height;
st_b = step_length;

st_a = sqrt(st_h^2 + (st_b/2)^2);
ang_st_A = atan(st_h/(st_b/2));
ang_st_B = pi - 2*ang_st_A;

%% Feasibility Check
if st_a >= leg_length
    x = [];
    return;
end

%% Leg Triangle
% If thigh length is equal to shank length, the defination is as same as
% the stance triangle. if not, the angle between the thigh an the stance
% triangle's left leg is defined as new angle 'A';

% base of leg triangle
leg_b = st_a;

% the vertex angle
ang_leg_B = acos((thigh_length^2 + shank_length^2 -leg_b^2) ...
    / (2*thigh_length*shank_length));

% The angle between the thigh and the stance triangle's left leg
ang_leg_A = acos((thigh_length^2 + leg_b^2 -shank_length^2) ...
    / (2*thigh_length*leg_b));

%% Position q
beta = ang_leg_B - pi;

% When standing with both feet together, the front foot is the left foot.
alpha_l = -phi - ang_st_B/2 + ang_leg_A;
alpha_r = -phi + ang_st_B/2 + ang_leg_A;

q = [phi; alpha_l; alpha_r; beta; beta];

%% Velocity q_dot
J = J_Hip_m(q, p);
dq = pinv(J) * [v_target; 0];

%% Assembly states
x = [q;dq];


end
