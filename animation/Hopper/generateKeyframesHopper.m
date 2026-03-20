function generateKeyframesHopper(data, epsilon, full3D, outputFolder, timeStep)
%GENERATEKEYFRAMESHOPPER Generates keyframes (snapshots) of a one-legged robot
%   Generates frames at fixed time intervals across all phases.
%   Final frame of each phase is saved separately.

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Determine 3D settings
    if full3D
        edges = false;
    else
        edges = true;
    end

    % Colors for different parts of the robot
    colors.torso     = [220, 220, 220]/255;
    colors.innerTube = [220, 220, 220]/255;
    colors.outerTube = [70, 70, 70]/255;
    colors.spring    = [10, 81, 131]/255;

    % Initialize the robot graphics
    hopper_graphic = OneLeggedHopper(full3D, edges, colors, epsilon);

    frameCounter = 1;
    t_global = 0; % global time reference

    % Loop over all phases
    for iDomain = 1:numel(data)
        tDomain = data(iDomain).t;
        tDomainShifted = tDomain + t_global;

        % Generate target times at fixed steps, **excluding the last time**
        t_targets = t_global:timeStep:(tDomainShifted(end) - 1e-12); 
        % small epsilon to avoid including the final frame

        for t_target = t_targets
            % Find closest index in this domain
            [~, idx] = min(abs(tDomainShifted - t_target));

            % Get state and input
            state_ = data(iDomain).x(idx, :)';
            input_ = data(iDomain).u(idx, :)';

            % Update robot graphics
            update(hopper_graphic, state_', input_', epsilon);

            % Save figure
            filename = sprintf('%s/keyframe_domain%d_frame%d.png', outputFolder, iDomain, frameCounter);
            saveas(gcf, filename, 'png');
            frameCounter = frameCounter + 1;
        end

        % --- Save the final frame of the phase separately ---
        state_ = data(iDomain).x(end, :)';
        input_ = data(iDomain).u(end, :)';
        update(hopper_graphic, state_', input_', epsilon);
        filename = sprintf('%s/keyframe_domain%d_final.png', outputFolder, iDomain);
        saveas(gcf, filename, 'png');

        % Update global time reference for next phase
        t_global = tDomainShifted(end) - timeStep; 
        % keep spacing consistent
    end

    fprintf('Keyframes saved to %s\n', outputFolder);
end
