function F_complete = getAnimationHopper5DOF(data,epsilon,full3D)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here

if full3D % change the 3D lighting effects
    edges         = false; % show edges for torso and legs
else
    edges         = true; % show edges for torso and legs
end
    
% Colors [R,G,B] for torso;leg;spring;
colors.torso     = [220, 220, 220]/255;
colors.innerTube = [220, 220, 220]/255;
colors.outerTube = [70, 70, 70]/255;
colors.spring    = [10,81,131]/255; 

% stateIn  = [x,y,phi,alpha,l]

hopper_graphic = OneLeggedHopper(full3D,edges,colors,epsilon);

F_complete = [];
for iDomain = 1:size(data,2)
    F(length(data(iDomain).t)) = struct('cdata',[],'colormap',[]); 
    frame_counter = 1;
    for i=1:length(data(iDomain).t)
        state_ = data(iDomain).x(i,:)';
        input_ = data(iDomain).u(i,:)';

        update(hopper_graphic, state_', input_', epsilon);
        % append
        F(frame_counter) = getframe(gcf);
        frame_counter = frame_counter + 1;
    end
    F_complete = [F_complete, F];

    clearvars F % s.t. F can be resized
end
end