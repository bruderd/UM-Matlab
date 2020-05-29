classdef Ksim
    %mpcsim: Class for performing simulations of fake systems
    %   Build from a system class (like arm) and an mpc class
    
    properties
        sys;
        mpc;
        scaledown;  % for scaling data into range [-1 , 1]
        scaleup;    % for unscaling data from range [-1 , 1]
        ref;        % reference trajectory
        shape;      % reference point or trajectory of shape parameters (for arm system only)
    end
    
    methods
        function obj = Ksim( system_class , mpc_class , varargin )
            %CLASS CONSTRUCTOR
            %   Detailed explanation goes here
            
            obj.sys = system_class; % hold onto entire class
            obj.mpc = mpc_class;    % hold onto entire class
            
            % copy the scaling functions for easier access
            obj.scaledown = mpc_class.scaledown;
            obj.scaleup = mpc_class.scaleup;
%             obj = obj.get_refscale;   % get the scale matrix for reference traj.
            
            % set default values of optional inputs
            ref = [];
            shape = [];
            
            % replace default values with user input values
            obj = obj.parse_args( varargin{:} );
        end
        
        % parse_args: Parses the Name, Value pairs in varargin
        function obj = parse_args( obj , varargin )
            %parse_args: Parses the Name, Value pairs in varargin of the
            % constructor, and assigns property values
            for idx = 1:2:length(varargin)
                obj.(varargin{idx}) = varargin{idx+1} ;
            end
        end
         
        %% Simulate the system with mpc controller
        
        % run_trial_mpc: Runs a simulation of system under linear mpc controller
        function results = run_trial_mpc( obj , ref , x0 , u0)
            %run_trial: Runs a simulation of system under mpc controller.
            %   Tries to follow the trajectory in ref and impose the
            %   shape constraints in shape_bounds.
            %   Assume ref and shape_bounds have same sampling frequency as
            %   sytem, and that they are already scaled to be consistent 
            %   with the lifted model.
            %   ref - struct containing reference trajectory with fields:
            %       t - vector of timesteps
            %       y - each row is a desired point at the corresponding timestep
            %   x0 - [1,nx] initial condtion
            %   u0 - [1,nu] initial input
            
            % shorthand
            nd = obj.mpc.params.nd;
            Np = obj.mpc.horizon;
            
            % set value of initial conditions to zero if none provided
            if nargin < 3
                x0 = zeros( nd+1 , obj.sys.params.nx );
                u0 = zeros( nd+1 , obj.sys.params.nu );
            elseif nargin < 4
                x0 = kron( ones( nd+1 , 1 ) , x0 );
                u0 = zeros( nd+1 , obj.sys.params.nu );
            else
                x0 = kron( ones( nd+1 , 1 ) , x0 );
                u0 = kron( ones( nd+1 , 1 ) , u0 );
            end
            y0 = obj.sys.get_y( x0 );
            
            % resample and scale the reference trajectory (TODO: FIX THIS AND THE REFERENCE TRAJECTORY FORMAT)
%             ref_Ts = obj.resample_ref( ref );
%             ref_sc = obj.scaledown.y( ref_Ts );
            ref_sc = obj.scaledown.ref( ref );
            
            % set initial condition
            initial.y = obj.scaledown.y(y0); initial.u = obj.scaledown.u(u0);
            [ initial , zeta0 ] = obj.mpc.get_zeta( initial );    % LINE NOT NEEDED
            
            % initialize results struct
            results = struct;
            results.T = [ 0 ];
            results.U = [ u0( end , : ) ];
            results.Y = [ y0( end , : ) ];
            results.K = [ 0 ];
            results.R = [ ref(1,:) ];
            results.X = [ x0( end , : ) ];
            results.Z = [ obj.mpc.lift.econ_full( zeta0' )' ]; % lifted states
            
            k = 1;
            while k < size( ref_sc , 1 )
                
                % current time
                t = k * obj.mpc.params.Ts;
                
                % get current state and input with delays
                if k == 1
                    current.y = obj.scaledown.y( y0 );   
                    current.u = obj.scaledown.u( u0 );  
                elseif k < nd + 1
                    y = [ y0( k : end-1 , : ) ; results.Y ];
                    u = [ u0( k : end-1 , : ) ; results.U ];
                    current.y = obj.scaledown.y( y );
                    current.u = obj.scaledown.u( u ); 
                else
                    y = results.Y( end - nd : end , : );
                    u = results.U( end - nd : end , : );
                    current.y = obj.scaledown.y( y ); 
                    current.u = obj.scaledown.u( u ); 
                end
                
                % isolate the reference trajectory over the horizon
                if k + Np <= size( ref_sc , 1 )
                    refhor = ref_sc( k : k + Np , :);
                else
                    refhor = ref_sc( k : end , : );     % repeat last entry
                end 
                
                % get optimal input over horizon
                if strcmp( obj.mpc.model_type , 'linear' )
                    [ U , z ] = obj.mpc.get_mpcInput( current , refhor );
                elseif strcmp( obj.mpc.model_type , 'bilinear' )
                    [ U , z ] = obj.mpc.get_mpcInput_bilinear( current , refhor );
                end
                
                % if a solution was not found, break out of while loop
                if any( isnan(U) )
                    break;
                end
                
                % isolate input for this step (may need to make it U(1,:)
                u_kp1_sc = U( 2 , : );
                
                % scaleup the input for the results
                u_kp1 = obj.scaleup.u( u_kp1_sc )';
                
%                 % simulate the system over one time-step (Using sysid model)
%                 z_k = z;
%                 u_k_sc = obj.scaledown.u( results.U(end,:) );  % need to use previously calculated input NEED TO MAKE THIS CLEANER!!!
%                 z_kp1 = obj.mpc.model.A * z_k + obj.mpc.model.B * u_k_sc';
%                 x_kp1 = obj.mpc.model.C * z_kp1;  % DON'T CARE ABOUT THIS JUST WANT IT TO WORK
%                 y_kp1_sc = x_kp1;  % output just scaled version of state since model was learned from observations
%                 y_kp1 = obj.scaleup.y( y_kp1_sc' )';  % scale output back up 
                
                % simulate the system over one time-step (using actual system model)
                x_k = results.X( end , : )';
                u_k = results.U(end,:)'; %current.u';
                x_kp1 = obj.sys.simulate_Ts( x_k , u_k , [0 0]' );
                y_kp1 = obj.sys.get_y( x_kp1' )';
                
                % record updated results
                results.T = [ results.T ; t ];
                results.U = [ results.U ; u_kp1' ];
                results.Y = [ results.Y ; y_kp1' ];
                results.K = [ results.K ; k ];
                results.R = [ results.R ; obj.scaleup.ref( ref_sc( k , : ) ) ];   % note that this is not scaled down
                results.X = [ results.X ; x_kp1' ];
                results.Z = [ results.Z ; z'  ]; % current lifted state
                
                k = k + 1;  % increment step counter
            end
        end
        

    end
end