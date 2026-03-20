function drawRABBIT(X, step, gamma, config)

import BF_floating.left_Stance.*;
% import BF_min.left_Stance.*;
% import BF_min_2.left_Stance.*;



persistent StanceFoot;
% StanceFoot = [0;0];
if isempty(StanceFoot)
    StanceFoot = [0;0];
end
% StanceFoot = [X(1);X(2)];
% stride = {1,2}
[p, ~] = getFullParameters(config);

% trafo = [cos(gamma) -sin(gamma); sin(gamma) cos(gamma)]';
% posStanceFoot = posStanceFootAUTO(x,p);
if step
    posStanceFoot = StanceFoot-posSwingFootAUTO(X,p);
    StanceFoot = posStanceFoot;
    % posStanceFoot = [0,0];
else
    posStanceFoot = StanceFoot;
end
% colors
grey = [0.2431 0.2667 0.2980];
brightblue = [0 0.7451 1];
darkblue   = [0 0.3176 0.6196];

% get coordinates of rigid bodies
posSwingFoot = posStanceFoot + posSwingFootAUTO(X,p);
posSwingKnee = posStanceFoot + posSwingKneeAUTO(X,p);
posStanceKnee = posStanceFoot + posStanceKneeAUTO(X,p);
posHip = posStanceFoot + posHipAUTO(X,p);
% disp(posHip);
posHead = posStanceFoot + posHeadAUTO(X,p);

% plot stance leg / tibia
plot([posStanceFoot(1) posStanceKnee(1)], [posStanceFoot(2) posStanceKnee(2)], '-', 'Color', darkblue, 'LineWidth', 6);
% plot stance leg / femur
plot([posStanceKnee(1) posHip(1)], [posStanceKnee(2) posHip(2)], '-', 'Color', darkblue, 'LineWidth', 6);
% % plot swing leg / femur
% plot([posHip(1) posSwingKnee(1)], [posHip(2) posSwingKnee(2)], '-', 'Color', darkblue, 'LineWidth', 6);
% % plot swing leg / tibia
% plot([posSwingKnee(1) posSwingFoot(1)], [posSwingKnee(2) posSwingFoot(2)], '-', 'Color', darkblue, 'LineWidth', 6);

% plot torso
plot([posHip(1) posHead(1)], [posHip(2) posHead(2)], '-', 'Color', grey, 'LineWidth', 6);

% plot swing leg / femur
plot([posHip(1) posSwingKnee(1)], [posHip(2) posSwingKnee(2)], '-', 'Color', brightblue, 'LineWidth', 6);
% plot swing leg / tibia
plot([posSwingKnee(1) posSwingFoot(1)], [posSwingKnee(2) posSwingFoot(2)], '-', 'Color', brightblue, 'LineWidth', 6);
% % plot stance leg / tibia
% plot([posStanceFoot(1) posStanceKnee(1)], [posStanceFoot(2) posStanceKnee(2)], '-', 'Color', brightblue, 'LineWidth', 6);
% % plot stance leg / femur
% plot([posStanceKnee(1) posHip(1)], [posStanceKnee(2) posHip(2)], '-', 'Color', brightblue, 'LineWidth', 6);
end