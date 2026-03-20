function getAnimationRABBIT(traj,gamma,config, filename)
%getAnimationRABBIT Summary of this function goes here
%   Detailed explanation goes here

X = traj.x;
ie = traj.ie;
snapShots = figure('Units', 'pixels', 'Position', [100, 100, 1200, 800]);
grid off; hold on; box off; axis off;
axis([-2, 2, -1, 2.1])
set(gcf, 'Color', 'w');  % Set figure background to white
set(gca, 'Color', 'w');  % Set axes background to white

% animation parameters
nframes = size(X,1); % # of frames to draw between [0, t] units of time

ax = gca; % autoscales by default
axis equal
grid off
set(ax,'XTickLabel',[]);
set(ax,'YTickLabel',[]);


% save frames throughout the animation
frames(nframes) = struct('cdata',[],'colormap',[]); 

% animation loop
idxFrame = 1;
j = 1;
% for stride = [1,2]
for i = 1:nframes
    cla(ax);  % Clear axes content before drawing
    % draw grounds
    plot(10*[-1 1], 10*tan(gamma)*[1 -1], 'k-', 'LineWidth', 3);
    % draw robot
    x = X(i,:)';
    if ~isempty(ie)
        if i == ie(j)
            j = j+1;
            step = true;
        else
            step = false;
        end
    else
        step = false;
    end
    drawRABBIT(x, step, gamma, config);
    axis([-2, 2, -1, 2.1]);
    % save drawing
    frames(idxFrame) = getframe(ax);
    idxFrame = idxFrame+1;
end
% end


if ~isempty(filename)
    movie2gif(frames, [filename,'.gif'],'DelayTime',.05,'LoopCount',Inf)
end
close(snapShots)
end
