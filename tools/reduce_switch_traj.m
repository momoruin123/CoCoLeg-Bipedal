function traj_min = reduce_switch_traj(traj)
% REDUCE_SWITCH_TRAJ Convert HF_noPitch trajectory to minimal coordinate representation
%
%   Switches phase order from [Stance, Flight] to [Flight, Stance] and
%   removes unnecessary coordinates for the minimal coordinate model.
%   Assumes input trajectory was generated using HF_noPitch model.
%
%   Inputs:
%     traj - trajectory structure with phases [Stance, Flight]
%
%   Outputs:
%     traj_min - modified trajectory with phases [Flight, Stance] and
%                reduced state dimensions for minimal coordinate model
%
%   author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    % Switch phase order: [Stance, Flight] -> [Flight, Stance]
    traj_min(1) = traj(2);
    traj_min(2) = traj(1);

    % Remove unnecessary coordinates from stance phase for minimal representation
    % Eliminates first 2 states (x, y positions) and states 5-6 (likely velocities)
    traj_min(2).x(:,1:2) = [];
    traj_min(2).x(:,5:6) = [];
end