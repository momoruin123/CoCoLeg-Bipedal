function drawRABBIT(config, X, phaseName, stride)
modelName = config.model_name;
LinkPositions_f = str2func([modelName, '.', phaseName,'.LinkPositions_f']);
[p, ~] = getFullParameters(config);

grey = [0.2431 0.2667 0.2980];
brightblue = [0 0.7451 1];
darkblue   = [0 0.3176 0.6196];

% Get links positions
links = LinkPositions_f(X,p);
links = links';
posStanceFoot = links(1,:);
posStanceKnee = links(2,:);
posHip        = links(3,:);
posHead       = links(4,:);
posSwingKnee  = links(5,:);
posSwingFoot  = links(6,:);

if stride == 1
    color_config = [darkblue; grey; brightblue];
    % plot stance leg / tibia
    plot([posStanceFoot(1) posStanceKnee(1)], [posStanceFoot(2) posStanceKnee(2)], '-', 'Color', color_config(1,:), 'LineWidth', 6);
    % plot stance leg / femur
    plot([posStanceKnee(1) posHip(1)], [posStanceKnee(2) posHip(2)], '-', 'Color', color_config(1,:), 'LineWidth', 6);
    % plot swing leg / femur
    plot([posHip(1) posSwingKnee(1)], [posHip(2) posSwingKnee(2)], '-', 'Color', color_config(3,:), 'LineWidth', 6);
    % plot swing leg / tibia
    plot([posSwingKnee(1) posSwingFoot(1)], [posSwingKnee(2) posSwingFoot(2)], '-', 'Color', color_config(3,:), 'LineWidth', 6);
    % plot torso
    plot([posHip(1) posHead(1)], [posHip(2) posHead(2)], '-', 'Color', color_config(2,:), 'LineWidth', 6);

elseif stride == 2
    color_config = [brightblue; grey; darkblue];
    % plot swing leg / femur
    plot([posHip(1) posSwingKnee(1)], [posHip(2) posSwingKnee(2)], '-', 'Color', color_config(3,:), 'LineWidth', 6);
    % plot swing leg / tibia
    plot([posSwingKnee(1) posSwingFoot(1)], [posSwingKnee(2) posSwingFoot(2)], '-', 'Color', color_config(3,:), 'LineWidth', 6);
    % plot stance leg / tibia
    plot([posStanceFoot(1) posStanceKnee(1)], [posStanceFoot(2) posStanceKnee(2)], '-', 'Color', color_config(1,:), 'LineWidth', 6);
    % plot stance leg / femur
    plot([posStanceKnee(1) posHip(1)], [posStanceKnee(2) posHip(2)], '-', 'Color', color_config(1,:), 'LineWidth', 6);
    % plot torso
    plot([posHip(1) posHead(1)], [posHip(2) posHead(2)], '-', 'Color', color_config(2,:), 'LineWidth', 6);

else
    error("Wrong stride idx!")
end

% Update axis
view_half_width = 1.5;
axis([posHip(1) - view_half_width, posHip(1) + view_half_width, -0.8, 2.2]);
% Display window
% axis([-2, 2, -1, 2.1]);
end
