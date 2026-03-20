function traj_full = trajToFulltraj(config, traj, isFullstride)
% Extend the traj, that only contains half of the stride, to the whole
% stride and floating base coordinates

% Initialization
for i = 1:numel(traj)   
    traj(i).idxStep = 1;
end
if isFullstride
    traj_full = repmat(traj, 1, 2);
else
    traj_full = traj;
end
nPhases = numel(traj_full);
isSecondStep_F = false;
isSecondStep_SS = false;

% Get model parameters
[p, ~] = getFullParameters(config);

% Extract configuration parameters
modelName  = config.model_name;
phaseNames = config.phaseSequence;
initPos_x0 = config.operatingCond.x0;   % Initial x-position [m]

for i = 1:nPhases
    traj_i = traj_full(i);
    currentPhase = traj_i.Phase;

    % Find the corresponding phase：Flight, SingleStance
    if strcmp(currentPhase, 'Flight')
        % Floating base coordinates are used in SingleStance phase
        X_f = traj_i.x;
       
        % Compute stance foot position
        % Get x-pos from the X_end of last phase
        X_end_f = traj_full(i-1).x(end, :);         % End states of last phases 
        x_end   = X_end_f(1);                       % x-pos of floating base
        X_f(:,1)  = X_f(:,1) + x_end;
        traj_i.x = X_f;
        
        % Swap legs
        if isSecondStep_F
            traj_i.idxStep = 2;
        end
        % Reverse tags from false to true
        isSecondStep_F = ~isSecondStep_F;

    elseif strcmp(currentPhase, 'SingleStance')
        % Minimal coordinates are used in SingleStance phase
        % Extract states from trajectory
        X_m  = traj_i.x;
        nq_m = size(X_m,2)/2;
        q_m  = X_m(:, 1:nq_m);
        dq_m = X_m(:, nq_m+1:end);
        % Initialization
        X_f = zeros(size(X_m,1), 14);

        % Compute stance foot position
        if i == 1
            % if it is the first phase, get x-pos from config
            extend_pos_x = initPos_x0;
        else
            % if not, get x-pos from the X_end of last phase
            pos_SwingFoot_f = str2func([modelName, '.', currentPhase, '.pos_SwingFoot_f']);
            X_end_f = traj_full(i-1).x(end, :);
            pos = pos_SwingFoot_f(X_end_f', p);
            extend_pos_x = pos(1);
        end

        % Extend to floating base trajectory
        for j = 1:size(X_m, 1)
            q_f  = [extend_pos_x, 0, q_m(j,:)];
            dq_f = [0, 0, dq_m(j,:)];
            X_f(j,:) = [q_f, dq_f];
        end
        traj_i.x = X_f;

        % Swap legs
        if isSecondStep_SS
            traj_i.idxStep = 2;
        end
        % Reverse tags from false to true
        isSecondStep_SS = ~isSecondStep_SS;
    else
        error("The trajectory contains an unsupported phase!")
    end

    traj_full(i) = traj_i;
end

end

function traj_swapped = swap_legs(traj_original)
% Swap legs to implement the rest half of the stride
% Extract data
X_swapped = traj_original.x;
u_swapped = traj_original.u;

% Swap states
X_swapped(:,[4,5])   = X_swapped(:,[5,4]);
X_swapped(:,[6,7])   = X_swapped(:,[7,6]);
X_swapped(:,[4,5]+7) = X_swapped(:,[5,4]+7);
X_swapped(:,[6,7]+7) = X_swapped(:,[7,6]+7);

% Swap input
u_swapped(:, [2,1]) = u_swapped(:, [1,2]);
u_swapped(:, [4,3]) = u_swapped(:, [3,4]);

% Summary
traj_swapped = traj_original;
traj_swapped.x = X_swapped;
traj_swapped.u = u_swapped;
end