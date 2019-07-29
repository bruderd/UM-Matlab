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
        
%         % dynamics as vector field
%         function  vf = 
        
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
            
        
    end
end











