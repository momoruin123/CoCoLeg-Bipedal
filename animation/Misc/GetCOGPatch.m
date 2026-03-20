% *************************************************************************
% function [f, v, c] = GetCOGPatch(x, y, phi, alpha, l, ualpha, ul, e)
%  
% Returns faces, vertices, and color-values representing a patch that shows
% the 3D view of the COG of a bipedal robot
%
% GetMonopod2DOFPatch(x, y, phi, alpha, ualpha, l, ul), without return
% arguments shows the patch in a new figure
%
% Input:  - x, y, phi (position and orientation of the main body)
%         - alpha, l (length and orientation of the leg)
%         - ualpha, ul (deflection of the actuators from their resting
%           position)
%         - e1, e2 (epsilon1 and epsilon 2 scaling factors)
%
% Output: - faces 'f', vertices 'v', and color values 'c' as used in the
%           'patch' function

function [f, v, c] = GetCOGPatch()
 
    %% COG patch as cylinder
    
    [x1,y1,z1] = cylinder([0,0.065,0.065,0],100); % cylinder with top and bottom surface
    z1([1,2],:)= 0; 
    z1([3,4],:)= 0.02;

    [f1,v1] = surf2patch(x1,y1,z1); % patch data of the cylinder from surface data
    c1 = ones(length(f1),3); % white colour
    % alternate quarters black colour
    c1(1:3:75,:) = repmat([0 0 0], 25, 1); 
    c1(151:3:225,:)=repmat([0 0 0], 25, 1); 
    c1(3:3:75,:) = repmat([0 0 0], 25, 1); 
    c1(153:3:225,:)=repmat([0 0 0], 25, 1);

    %rotate 90 around x axis
    v1 = TransformVertices(v1, [1 0 0;0 0 -1; 0 1 0], [0 0 0]);
 
    % black border for COG
    [x2,y2,z2] = cylinder([0,0.07,0.07,0],100); % cylinder with top and bottom surface
    z2([1,2],:)= 0; 
    z2([3,4],:)= 0.015;
    [f2,v2] = surf2patch(x2,y2,z2); % patch data of the cylinder from surface data
    c2 = zeros(length(f2),3); % white colour
    v2 = TransformVertices(v2, [1 0 0;0 0 -1; 0 1 0], [0 -0.0025 0]);
    [v,f,c]=AddPatchesWithColor(v1,f1,c1,v2,f2,c2); 

%     %% Epsilon values for COG positions
%     
%     % Legs
%     dl = epsilon1*0.3; % COG position for legs
%     % Torso
%     if epsilon1 == 0
%         dt = 0;
%     else
%         dt =  (epsilon1*0.35); % COG position for torso which has min height of (2*0.2)
%     end
%     theta_t = epsilon2+0.4; % COG scaling of torso according to epsilon 2 values
%     legScale = 0.05+epsilon2*0.35; % COG scaling of torso according to epsilon 2 values
%     
%     %% Main body COG
%     v_ = TransformVertices(v1,Body312dc([0,0,-phi]),[0,0,0]); % rotate
%     v_ = TransformVertices(v_,diag([theta_t,1,theta_t]),[0,0,0]); % Inertia scaling
%     v_ = TransformVertices(v_,eye(3),[dt*sin(-phi),0,dt*cos(-phi)]); % move with torso
%     v = TransformVertices(v_, eye(3), [x,-0.1,y]); % shift, COG infront of leg
%     f = f1;
%     c = c1;
%         
%     %% Left leg COG --- (left leg is in the back)
%     v2 = TransformVertices(v1,diag([legScale,1,legScale]),[0,0,0]); % Smaller COG marker for legs  
%     
%     v_ = TransformVertices(v2,Body312dc([0,0,-alpha_L-phi]),[0,0,0]); % rotate
%     v_ = TransformVertices(v_,eye(3),[dl*sin(alpha_L+phi),0,-dl*cos(alpha_L+phi)]); % move with leg
%     vL2 = TransformVertices(v_,eye(3), [x,+0.05,y]); % shift
%     [v,f,c]=AddPatchesWithColor(v,f,c,vL2,f1,c1);
%     
%     %% Right leg COG --- (right leg is in the front)   
%     v_ = TransformVertices(v2,Body312dc([0,0,-alpha_R-phi]),[0,0,0]); % rotate
%     v_ = TransformVertices(v_,eye(3),[dl*sin(alpha_R+phi),0,-dl*cos(alpha_R+phi)]); % move with leg
%     vR2 = TransformVertices(v_,eye(3), [x,-0.3,y]); % shift
%     [v,f,c]=AddPatchesWithColor(v,f,c,vR2,f1,c1); 
    
end
        
        
        
        