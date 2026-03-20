% *************************************************************************
%
% classdef OutputCLASS
%
% This abstract class defines the interface used for the (graphical) output
% of the simulations.  It defines functions for initialization and regular
% updating.  Specific implementations could reach from simple line graphs
% to more complex 2D and 3D representations. It should also be possible to
% plot the states as a function of time, or save them to disc for later
% use.
%
% Properties: - 'slowDown' Determines the speed of output. '1' is real time,
%                          '0' as fast as possible.
%             - 'rate'     Determines the update rate (in units of
%                          simulation time)
% Methods:    - t = getTimeVector(obj, tStart, tEnd) Tells the calling
%                          function for which point in time output is
%                          required. 
%             - obj = update(obj, y, z, t, u) Called to update the (screen)
%                          output
%
%   C. David Remy remy@inm.uni-stuttgart.de
%   Matlab R2018
%   4/22/2020
%   v21
%
% Based on the paper:
%  'A MATLAB Framework For Gait Creation', 2011, C. David Remy, Keith
%  Buffinton, and Roland Siegwart,  International Conference on
%  Intelligent Robots and Systems, September 25-30, San Francisco, USA 

classdef OutputCLASS
    properties 
        slowDown; 
        % This value determines the rate of the output. 1 is realtime.
        % bigger values slow it down according to the given factor, lower 
        % values speed it up. 0 Means as fast as possible.  
        rate;
        % This value (in units of time) determines the time rate at which 
        % output is created. 
    end
    methods
        % standard constructor (this should be overwritten to allow the
        % initialization of the graphic objects:
        function obj = OutputCLASS()
            obj.slowDown = 1;     % run in real-time 
            obj.rate     = 0.04;  % with 25 fps
        end
        % This function is used to pass the desired refresh rate to the
        % simulation:
        % Input:  - The beginning of the simulation 'tStart'
        %         - The end of the simulation 'tEnd'
        % Output: - The time grid vector 't' on which output is produced.
        %           It will be passed on to the integrator, which in
        %           turn calls the 'update' function.    
        function t = getTimeVector(obj, tStart, tEnd)
            tEnd = min(tEnd, 1e3); % limit the time, as infinite integration is not possible
            t = [tStart:obj.rate:tEnd,tEnd];
            t = [t(diff(t)~=0), t(end)];
        end
    end
    methods (Abstract)
        % This function is called every time the data has changed and
        % needs to be displayed.  It needs to be rewritten for every
        % specific implementation.
        % Input:  - The continuous state vector y
        %         - The discrete state vector z
        %         - The current time t
        %         - The control state vector u (u = [] for passive systems)
        obj = update(obj, y, z, t, u, x_st)
    end
end
% *************************************************************************
% *************************************************************************