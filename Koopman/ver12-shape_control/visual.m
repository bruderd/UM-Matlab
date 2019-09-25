classdef visual
    %visual: This class stores functions used for visualizing 2D/3D polys.
    %   
    
    properties
        Property1
    end
    
    methods
        function obj = visual(varargin)
            %visual: Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
        end
    
        % parse_args: Parses the Name, Value pairs in varargin
        function obj = parse_args( obj , varargin )
            %parse_args: Parses the Name, Value pairs in varargin of the
            % constructor, and assigns property values
            for idx = 1:2:length(varargin)
                obj.(varargin{idx}) = varargin{idx+1} ;
            end
        end
        
    end

    methods(Static)    
       
        % plot_fromCoeffs: Plots 2D/3D polynomials specified as coefficients
        function [ shape , pshape , fig ] =  plot_fromCoeffs( coeffs , fig )
            % plot_fromCoeffs: Plots 2D/3D polynomials specified as coefficients
            %   coeffs - [ dim x degree ] matrix. dim is the dimension of 
            %     the space (2D or 3D), degree is the degree of the polys
            %     Note that coeffs are written using the following
            %     convention:
            %       [ a b c d ... ] --> ax + bx^2 + cx^3 + ...
         
            if nargin < 2
                fig = figure;
                % set axes and other preferences ...
            end
            
            % if 2D coeffs are provided, add zero y-coefficients
            if( size(coeffs,1) == 2 ) 
                coeffs = [ coeffs(1,:) ;...
                           zeros(1,size(coeffs,2)) ;...
                           coeffs(2,:) ];
%                 view(2); % change view of the plot
            end
            
            % get the shape
            px = fliplr( [ 0 , coeffs(1,:) ] );
            py = fliplr( [ 0 , coeffs(2,:) ] );
            pz = fliplr( [ 0 , coeffs(3,:) ] );
            
            pol_x = zeros(101,1); pol_y = zeros(101,1); pol_z = zeros(101,1);
            for i = 1:101
                pol_x(i) = polyval(px,0.01*(i-1));
                pol_y(i) = polyval(py,0.01*(i-1));
                pol_z(i) = polyval(pz,0.01*(i-1));
            end
            shape = [ pol_x , pol_y , pol_z];
            
            % plot the shape
            pshape = plot3( pol_x , pol_y , pol_z , 'b' , 'LineWidth' , 3);  
            
            % plot the base of the robot
            hold on;
            patch( [-1 1 1 -1] , [ -1 -1 1 1 ] , 'k' , 'FaceAlpha' , 0.1)
            hold off;
            
            % change the viewpoint
            grid on; box on;
            fig.CurrentAxes.ZDir = 'Reverse';
%             fig.CurrentAxes.CameraPosition = [-9 -12 +7];
            view(-45,15)
        end
        
        % plot_fromMarkers: Plots 2D/3D polynomials specified as a list of marker positions
        function [ shape , pshape , pmarkers , fig ] = plot_fromMarkers( degree , points , positions , orient , fig )
            % plot_fromMarkers: Plots 2D/3D polynomials specified as a list of marker positions
            
            % set optional argument default values
            if nargin < 4
                orient = [];    
                fig = figure;
                % set axes and other preferences ...
            elseif nargin < 5
                fig = figure;
                % set axes and other preferences ...
            end
            
            % get shape coefficients
            [ coeffs , ~ ] = visual.points2poly(degree, points, positions, orient);
            
            % plot the shape
            [ shape , pshape , fig ] =  visual.plot_fromCoeffs( coeffs , fig );
            
            % plot the marker points
            if( size(points,2) == 2 ) 
                points = [ points(:,1) ,...
                           zeros(size(points,1),1) ,...
                           points(:,2) ];
%                 view(0,0);
            end
            hold on;
            pmarkers = plot3( points(:,1) , points(:,2) , points(:,3) , 'r.' , 'MarkerSize' , 30);
            hold off;
        end
        
        % points2poly (convert marker points into a polynomial)
        function [coeffs, obs_matrix] = points2poly(degree, points, positions, orient)
            %points2poly: Finds polynomial that best goes through a set of points.
            % Polynomial starts at the origin, and its domain is [0,1].
            % "Resting configuration" is along the yaxis (2d) or zaxis (3d)
            %   degree - scalar, maximum degree of the polynomial
            %   points - matrix, rows are xy(z) coordinates of points
            %   positions - (row) vector of [0,1] positions of points on the arm.
            %   orient - (row) vector, orientation of the the end effector (complex number for 2d case, quaternion for 3d case)
            %   coeffs - matrix, rows are the coefficients of the polynomial in each
            %            coordinate where [a b c ...] -> ax^1 + bx^2 + cx^2 + ...
            %   obs_matrix - matrix that converts from state vector to coefficients
            
            % set default for optional argument(s)
            if nargin < 4
               orient = [];
            end
            
%             % use this for checking orient when its provided
%             if all( size(orient) ~= [1,2] )
%                 error('orientation for 2d system must be a complex number specified as [1x2] vector');
%             end
            
            % for the 2d case (will consider the 3d case later)
            if size(points,2) == 2

                % generate virtual points to provide slope constraint at the base & end
                startpoint = [ 0 , 1e-2 ];
                if isempty( orient )
                    positions_supp = [ 0 , 1e-2 , positions];
                    points_supp = [0 , 0 ; startpoint ; points];
                else
                    positions_supp = [ 0 , 1e-2 , positions , 1+1e-2 ];
                    endpoint = ( orient/norm(orient) )*1e-2 + points(end,:);    % add point extended from last link
                    points_supp = [0 , 0 ; startpoint ; points ; endpoint];
                end
                
                
                % generate A matrix for least squares problem
                A = zeros( length(positions_supp) , degree );
                for i = 1 : degree
                    A(:,i) = reshape( positions_supp , [] , 1) .^ i;
                end
                
                % separate x and y corrdinates of points
                points_x = points_supp(:,1);
                points_y = points_supp(:,2);
                
                % find polynomial coefficients
                obs_matrix = pinv(A);
                coeffs_vec_x = obs_matrix * points_x;
                coeffs_vec_y = obs_matrix * points_y;
                
                % make coeffs a matrix to where each row is coeffs for one dimension
                coeffs = [ coeffs_vec_x' ; coeffs_vec_y' ];
            else
                error('points matrix must be nx2');
            end
        end
        
        % animate
        function animate( t , varargin )
            %animate: a sequence of 3D polynomials in time
            %   t - time vector from simulation
            %   varargin = Name, Value pairs
            %       Names:
            %           -'markers' 
            %           -'coeffs'
            %           -'name'
            %           -'degree'
            %       Values:
            %           -struct with the following fields
            %               points - [length(t) x 3degree]
            %               positions - [1 x numMarkers] location of markers in parametrized coordinates 
            %               orient - [lenght(t) x 4] orientation of end effector as quaternion (eventually)
            %           [ length(t) x 3degree ] - polynomial coeffs
            %           [ string ] - name of video file
            
            input.name = 'test-3Danimation';    % default filename
            input.degree = 3;   % default polynomial degree
            
            % sort out the optional inputs
            for idx = 1:2:length(varargin)
                input.(varargin{idx}) = varargin{idx+1} ;
            end
            
            fig = figure;   % create figure for the animation
            axis([-1 , 1 , -1 , 1 , -1 , 1] * 0.5)
            set(gca,'Zdir','reverse')
            xlabel('x(m)')
            ylabel('y(m)')
            zlabel('z(m)')
            
            % Prepare the new file.
            vidObj = VideoWriter( ['animations' , filesep , input.name , '.mp4'] , 'MPEG-4' );
            vidObj.FrameRate = 30;
            open(vidObj);
            
            set(gca,'nextplot','replacechildren', 'FontUnits' , 'normalized');
            
            totTime = t(end);    % total time for animation (s)
            nsteps = length(t); % total steps in the simulation
            totFrames = 30 * totTime;   % total frames in 30 fps video
            
            % run animation fram by frame
            for i = 1:totFrames
                
                index = floor( (i-1) * (nsteps / totFrames) ) + 1;   % skips points between frames
                
                hold on;
                if isfield( input , 'markers' ) % get shape from markers
                    points_vec = input.markers.points(index,:);
                    points_mat = reshape( points_vec' , [ length(points_vec) / length(input.markers.positions) , length(input.markers.positions) ] )';
                    [ ~ , p1 , m1 , ~ ] = visual.plot_fromMarkers( input.degree , points_mat , ...
                                      input.markers.positions, input.markers.orient , fig );    % need to fix orient argument
                end
                if isfield( input , 'coeffs' ) % get shape from coeffs
                    coeffs_mat = [ coeffs.x(index,:) ; coeffs.y(index,:) ; coeffs.z(index,:) ];
                    [ ~ , p2 , ~ ] = visual.plot_fromCoeffs( coeffs_mat , fig );
                end
                hold off;
                
                % write each frame to the file
                currFrame = getframe(fig);
                writeVideo(vidObj,currFrame);
                
                % clear the plots for next frame
                if exist('p1','var')
                    delete(p1);
                    delete(m1);
                end
                if exist('p2','var')
                    delete(p2);
                end
            end
            close(vidObj); 
        end
        
    end
end

















