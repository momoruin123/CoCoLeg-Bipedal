function getAnimationRABBIT(traj,gamma,config, filename)

% Extract traj for frames
nframes = 0; % # of frames to draw between [0, t] units of time
for i = 1:size(traj,2)
    nframes = nframes + size(traj(i).x, 1);
end

snapShots = figure('Units', 'pixels', 'Position', [100, 100, 1200, 800]);
grid off; hold on; box off; axis off;
axis([-2, 2, -1, 2.1])
set(gcf, 'Color', 'w');  % Set figure background to white
set(gca, 'Color', 'w');  % Set axes background to white

ax = gca; % autoscales by default
axis equal
grid off
set(ax,'XTickLabel',[]);
set(ax,'YTickLabel',[]);


% save frames throughout the animation
frames(nframes) = struct('cdata',[],'colormap',[]); 

% animation loop
idxFrame = 1;

for phase = 1:size(traj,2)
    x = traj(phase).x;
    phaseName = traj(phase).Phase;
    idxStep = traj(phase).idxStep;
    for i = 1:size(x,1)
        cla(ax);  % Clear axes content before drawing
        % draw grounds
        plot(10*[-1 1], 10*tan(gamma)*[1 -1], 'k-', 'LineWidth', 3);
        % draw robot
        x_i = x(i,:)';
        drawRABBIT(config, x_i, phaseName, idxStep);
        % save drawing
        frames(idxFrame) = getframe(ax);
        idxFrame = idxFrame+1;
    end
end
% end



if ~isempty(filename)
    movie2gif(frames, [filename,'.gif'],'DelayTime',.05,'LoopCount',Inf)
end
close(snapShots)
end
