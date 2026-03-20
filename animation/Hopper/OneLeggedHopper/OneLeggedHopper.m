% *************************************************************************
%
% classdef PrismaticMonopod3DCLASS(full3D) < OutputCLASS
%
% Three dimensional graphics of a prismatic monopod with series elastic
% actuation.  In the graphical representation, the robot is drawn with two
% legs which move together.  Since the motion is restricted to be planar,
% this is identical to a monopod.
%
% The system is initialized with a Boolean flag that indicates whether the
% view is fully three dimensional, or rather flat from the side.  The
% motion of the model is, however, in both cases 2D. 
%
%
% Properties: - NONE
% Methods:    - NONE
%
%
% Created by C. David Remy on 03/14/2011
% MATLAB 2010a
%
% Documentation:
%  'A MATLAB Framework For Gait Creation', 2011, C. David Remy (1), Keith
%  Buffinton (2), and Roland Siegwart (1),  International Conference on
%  Intelligent Robots and Systems, September 25-30, San Francisco, USA 
%
% (1) Autonomous Systems Lab, Institute of Robotics and Intelligent Systems, 
%     Swiss Federal Institute of Technology (ETHZ) 
%     Tannenstr. 3 / CLA-E-32.1
%     8092 Zurich, Switzerland  
%     cremy@ethz.ch; rsiegwart@ethz.ch
%
% (2) Department of Mechanical Engineering, 
%     Bucknell University
%     701 Moore Avenue
%     Lewisburg, PA-17837, USA
%     buffintk@bucknell.edu
%
%   See also OUTPUTCLASS.
%
classdef OneLeggedHopper < OutputCLASS 
    properties (SetAccess = 'private', GetAccess = 'private')
        fig; % The output window
        % Patch objects used in the 3D representation
        HopperPatch;
        FloorPatch;
        % Vertices to which the transformation will be applied
        FloorVertices;
        % homotopy parameter
        epsilon
        % View Flag
        full3D
        edges
        % patch colors
        colors
    end
    methods
        function obj = OneLeggedHopper(full3D,edges,colors,epsilon)
            obj.slowDown = 1; % Run this in real time.
            obj.rate     = 0.04;   % with 25 fps
            obj.full3D   = full3D;
            obj.epsilon  = epsilon;
            obj.colors   = colors;
            obj.edges    = edges;
            % Initialize the 3D graphics
            obj.fig = figure();
            clf(obj.fig);
            % Set some window properties
            set(obj.fig,'Name','Hopper Animation');  % Window title
            set(obj.fig,'Color','w');         % Background color
            
            %% Create the Objects:  
            % Create the floor:
            [f,v,c] = Get2DGroundPatch();
            p = patch('faces', f, 'vertices', v, 'FaceVertexCData', c, 'FaceColor', 'flat');
           
            set(p, 'FaceLighting','none'); % Set the renderer none, so pure white colour is achieved
            set(p, 'FaceColor','flat');     % Set the face coloring
            set(p, 'EdgeColor','none');     % Don't show edges
            % Store theses objects for later use (the ground must be
            % shifted later, if the hopper moves)
            obj.FloorPatch    = p;
            obj.FloorVertices = v;
            
            % Create the hopper (with some arbitrary initial condition)
            x     = 0; 
            y     = 1.2;
            phi   = 0.1;
            alpha = -0.1;
            l     = 1;
            u     = -[0.4;0.2];
            [f,v,c] = GetOneLeggedHopperPatch(x,y,phi,alpha,l,u,obj.epsilon,obj.edges,obj.colors);
            p = patch('faces', f, 'vertices', v, 'FaceVertexCData', c, 'FaceColor', 'flat');
            if full3D
                set(p, 'FaceLighting','phong');  % Set the renderer
            else
                set(p, 'FaceLighting','none');  % Set the renderer
            end
            set(p, 'FaceColor','flat');      % Set the face coloring
            set(p, 'EdgeColor','none');      % Don't show edges
            obj.HopperPatch = p;
            %if COG
%                 [f,v,c] = GetCOGPatch();
%                 p = patch('faces', f, 'vertices', v, 'FaceVertexCData', c, 'FaceColor', 'flat');
%                 set(p, 'FaceLighting','none');  % Set the renderer
%                 set(p, 'FaceColor','flat');      % Set the face coloring
%                 set(p, 'EdgeColor','none');      % Don't show edges
%                 %obj.COGPatch = p;
            %end
            %% set up view:
            axis off
            box off
            axis equal
            camproj('orthographic');
%                 camtarget([0, 0, +1.2])
%                 campos ([0, -10, +1.7]); 
            camtarget([0, 0, +1.1])
            campos ([0, -7.9, +1.1]);
            camup ([0,0,1]);
            camva(15)   %~ updated camva(10)
            
            %% Create illumination:
            light('Position',[5 -15 12],'Style','local');   % parallel light, coming in the direction given in 'Position'
        end
        function obj = update(obj, state, tau, epsilon)
            % rewrite state in more usable names
            x     = state(1);
            y     = state(2);
            phi   = state(3);
            alpha = state(4);
            l     = state(5);
            u     = tau(:);

            % The floor is shifted to multiples of its pattern, so it's
            % not noticeable in the final graphics: 
            v = TransformVertices(obj.FloorVertices,...
                                  diag([1,1,1]),...
                                  [floor(x/2)*2,0,-0.05]);
            set(obj.FloorPatch,'Vertices',v); % Main body    
            % The hopper:
            [~, v] = GetOneLeggedHopperPatch(x,y,phi,alpha,l,u,epsilon,obj.edges,obj.colors);
            set(obj.HopperPatch,'Vertices',v); 
            % Set camera:
            camtarget([x, 0, +1.1])
            campos ([x, -7.9, +1.1]);
            drawnow();
        end
    end
end