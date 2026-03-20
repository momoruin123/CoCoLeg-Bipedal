function traj = flow(config, dt)
% FLOW wrapper function for the flow functions specific for every trajectory.
%   traj = flow(config, dt) calls the flow function to run a simulation the system dynamics through the 
%   specified phase sequence using the given configuration and time step
%
%   Inputs:
%     config - system configuration structure
%     dt     - time step or vector of time steps for each phase
%
%   Output:
%     traj   - trajectory structure containing states, controls, time, 
%              parameters, and constraint forces

     if startsWith(config.model_name, 'h', 'IgnoreCase', true)
         traj = flow_hopper(config, dt);
     elseif startsWith(config.model_name, 'b', 'IgnoreCase', true)
         traj = flow_biped(config, dt);
     else
        error('Class not supported yet. Please adapt code to accomodate your class');
     end
end