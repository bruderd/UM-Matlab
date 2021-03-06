classdef sys
    %sys: Class with system properties and equations of motion inside
    %   This is for a planar manipulator arm specifically.
    
    properties
        params struct;  % model parameters
        fcns struct;   % equation of motion functions
    end
    
    methods
        function obj = sys( params , fcns )
            %Construct an instance of EOM class
            %   Detailed explanation goes here
            obj.params = params;
            obj.fcns = fcns;
        end
        
        %% transformations
        % alpha2theta
        function theta = alpha2theta( obj , alpha )
            theta = alpha2theta( alpha );
        end
        
        % alpha2xvec (gives x vectorized)
        function [ x , xcm ] = alpha2xvec( obj , alpha )
            [ x ,xcm ] = alpha2x( alpha , obj.params);
        end
        
        % alpha2x (gives x where rows are x,y coordinate pairs)
        function [ x , xcm ] = alpha2x( obj , alpha )
            [ x_vec ,xcm_vec ] = alpha2x( alpha , obj.params);
            x = reshape( x_vec , [ 2 , obj.params.Nlinks+1 ] )';
            xcm = reshape( xcm_vec , [ 2 , obj.params.Nlinks ] )';
        end
        
        
        %% equations of motion
        
        % get_massMatrix
        function Dq = get_massMatrix( obj , alpha )
            Dq = obj.fcns.get_massMatrix(alpha);
        end
        
        % get_nonInert
        function nonInert = get_nonInert( obj , alpha , alphadot , u )
            nonInert = obj.fcns.get_nonInert( alpha , alphadot , u );
        end
        
        % get_damp
        function damp = get_damp( obj , alphadot )
            damp = obj.fcns.get_damp( alphadot );
        end
        
        % get_input
        function input = get_input( obj , alpha , u )
            input = obj.fcns.get_damp( alpha , u );
        end
        
        % get_dampNinput
        function danpNinput = get_dampNinput( obj , alpha , alphadot , u )
            dampNinput = obj.fcns.get_dampNinput( alpha , alphadot , u );
        end
        
        % vf (dynamics as vector field)
        function Alphadot = vf( obj , t , Alpha , u )
            %vf: Explicit dynamics for robot
            %   Alpha = [ alpha ; alphadot ];
            %   Alphadot = [ alphadot ; alphaddot ];
            
            params = obj.params;
            
            alpha = Alpha( 1 : params.Nlinks );
            alphadot = Alpha( params.Nlinks+1 : end );
            
            Dq = obj.get_massMatrix( alpha );
            nonInert = obj.get_nonInert( alpha , alphadot , u );
            
            % solve for acceleration terms
            alphaddot = - Dq \ nonInert;
            
            % define output
            Alphadot = [ alphadot ; alphaddot ];
        end
        
        %% sensing
        
        % get_markers (simulated mocap)
        function markers = get_markers( obj , alpha )
            [ x , ~ ] = obj.alpha2x( alpha );
            markers = x( 1 : obj.params.nlinks : end , : );
        end
        
        % get_shape
        function [ shape , coeffs ] = get_shape( obj , alpha , degree)
            points = get_markers( obj , alpha );   % coordinates of mocap markers
            positions = obj.params.markerPos;    % relative location of markers on arm [0,1]
            theta = alpha2theta( alpha );
            orient = theta2complex( theta );    % orientaton of end effector
            coeffs = points2poly( degree , points , positions , orient );    % convert points of a polynomial
            
            % get the shape
            px = fliplr( [ 0 , coeffs(1,:) ] );
            py = fliplr( [ 0 , coeffs(2,:) ] );
            
            pol_x = zeros(101,1); pol_y = zeros(101,1);
            for i = 1:101
                pol_x(i) = polyval(px,0.01*(i-1));
                pol_y(i) = polyval(py,0.01*(i-1));
            end
            shape = [ pol_x , pol_y ];
        end
        
        %% visualization
        
        % def_fig (defines a default figure for plotting arm
        function fig = def_fig( obj )
            % set up figure
            fig = figure;
            axis([-obj.params.L, obj.params.L, -0.5*obj.params.L, 1.5*obj.params.L])
            set(gca,'Ydir','reverse')
            xlabel('x(m)')
            ylabel('y(m)')
        end
        
        % plot_arm
        function ph = plot_arm( obj , alpha )
            % convert to xy-coordinates
            [ X , ~ ] = alpha2x(alpha, obj.params);
            x = [0; X(1:2:end)];
            y = [0; X(2:2:end)];
            
            % add markers
            markers = obj.get_markers( alpha );
            
            % plot it
            hold on
            ph(1) = plot( x, y, '-o' );
            ph(2) = plot( markers(:,1) , markers(:,2) , 'r*');
            hold off
        end
        
        % plot_arm_shape
        function ph = plot_arm_shape( obj , alpha , degree )
            % plot the arm
            ph = obj.plot_arm( alpha );
            
            % get the shape
            [ shape , ~ ] = obj.get_shape( alpha , degree );
            
            % plot it
            hold on
            ph(3) = plot( shape(:,1) , shape(:,2) , 'r');
            hold off
        end
        
        % animate_arm
        function animate_arm( obj, t , y , varargin)
            %animate_arm: Animate a simualtion of the arm
            %   t - time vector from simulation
            %   y - state vector from simulation (alpha and alphadot)
            %   varargin{1} = degree - degree of the shape polynomial (default: 3)
            %   varargin{2} = name - name of the video file (default: sysName)
            
            % deal with optional inputs
            if length(varargin) == 2
                degree = varargin{1};
                name = varargin{2};
            elseif length(varargin) == 1
                degree = varargin{1};
                name = obj.params.sysName;
            else
                degree = 3;
                name = obj.params.sysName;
            end
            
            alpha = y(: , 1:obj.params.Nlinks );   % joint angles over time
            
            fig = figure;   % create figure for the animation
            axis([-obj.params.L, obj.params.L, -0.5*obj.params.L, 1.5*obj.params.L])
            set(gca,'Ydir','reverse')
            xlabel('x(m)')
            ylabel('y(m)')
            
            % Prepare the new file.
            vidObj = VideoWriter( ['animations' , filesep , name , '.mp4'] , 'MPEG-4' );
            vidObj.FrameRate = 30;
            open(vidObj);
            
            set(gca,'nextplot','replacechildren', 'FontUnits' , 'normalized');
            
            totTime = t(end);    % total time for animation (s)
            nsteps = length(t); % total steps in the simulation
            totFrames = 30 * totTime;   % total frames in 30 fps video
            
            % run animation fram by frame
            for i = 1:totFrames
                
                index = (i-1) * floor( nsteps / totFrames ) + 1;   % skips points between frames
                
                [ X , ~ ] = obj.alpha2x( alpha(index,:)' );
                x = X(:,1);
                y = X(:,2);
                marker = obj.get_markers( alpha(index,:) );   % get mocap sensor location
                [shape , ~ ] = obj.get_shape( alpha(index,:) , degree); % get polynomial approx of shape (3rd order)
                
                hold on;
                p1 = plot(x, y, 'b-o');
                p2 = plot( marker(:,1) , marker(:,2) , 'r*');
                p3 = plot( shape(:,1) , shape(:,2) , 'r');
                hold off;
                
                % write each frame to the file
                currFrame = getframe(fig);
                writeVideo(vidObj,currFrame);
                
                delete(p1); delete(p2); delete(p3);
            end
            
            close(vidObj);
        end
            
        %% simulation
        
        % simulate system under random "ramp and hold" inputs
        function [ t , y , u ] = simulate_arm( obj , tf , Tramp )
            %simulate: simulate system under random "ramp and hold" inputs
            %   tf - length of simulation(s)
            %   Tramp - ramp period length
            
            % time steps
            t = ( 0 : obj.params.Ts : tf )';    % all timesteps
            tswitch = ( 0 : Tramp : tf )';  % input switching points
            
            % table of inputs
            numPeriods = ceil( length(tswitch) / 2 );
            inputs_nohold = obj.params.umax .* ( 2*rand( numPeriods , obj.params.Nmods ) - 1 );  % table of random inputs
            inputs_hold = reshape([inputs_nohold(:) inputs_nohold(:)]',2*size(inputs_nohold,1), []); % repeats rows of inputs so that we get a hold between ramps
            u = interp1( tswitch , inputs_hold( 1:length(tswitch) , : ) , t );
            
            % initial condition (resting)
            a0 = zeros( obj.params.Nlinks , 1 );
            adot0 = zeros( obj.params.Nlinks , 1 );
            
            % simulate system
            [ y ] = ode5( @(t,x) obj.vf( t , x , u( floor(t/obj.params.Ts) + 1 , : )' ) , t , [ a0 ; adot0 ] );  % with numerical inversion, fixed time step
        end
           
        
    end
end











