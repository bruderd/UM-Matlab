classdef mpcsim
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
        function obj = mpcsim( system_class , mpc_class , varargin )
            %CLASS CONSTRUCTOR
            %   Detailed explanation goes here
            
            obj.sys = system_class; % hold onto entire class
            obj.mpc = mpc_class;    % hold onto entire class
            
            % copy the scaling matrices for easier access
            obj.scaledown = mpc_class.params.scaledown;
            obj.scaleup = mpc_class.params.scaleup;
            obj = obj.get_refscale;   % get the scale matrix for reference traj.
            
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
        
        %% Define reference trajectories and state constraints
        
        % get_refscale: creates matrices for scaling the reference trajectory
        function obj = get_refscale( obj )
            % get_refscale: creates matrices for scaling the reference trajectory
            
            scaledown_vec = obj.mpc.projmtx * diag( obj.scaledown.z );     % gets scalefactor for each element of the reference trajectory
            scaleup_vec = obj.mpc.projmtx * diag( obj.scaleup.z );     % gets scalefactor for each element of the reference trajectory
            
            obj.scaledown.ref = diag( scaledown_vec );
            obj.scaleup.ref = diag( scaleup_vec );
        end
        
        %% Simulate a the system with mpc controller
        
        % run_trial: Runs a simulation of system under mpc controller
        function results = run_trial( obj , ref , shape_bounds )
            %run_trial: Runs a simulation of system under mpc controller.
            %   Tries to follow the trajectory in ref and impose the
            %   shape constraints in shape_bounds.
            %   Assume ref has same sampling frequency as sytem.
            
            % shorthand
            nd = obj.mpc.params.nd;
            Np = obj.mpc.horizon;
            
            % scale the reference trajectory
            ref_sc = ref * obj.scaledown.ref;
%             ref_sc = ref;   % reference is already scaled
            
            % set initial condition (may add an input argument for this later)
            x0 = zeros( nd+1 , obj.sys.params.nx );
            y0 = obj.sys.get_y( x0 );   % get corresponding outputs
            u0 = zeros( obj.mpc.params.nd+1 , obj.sys.params.nu );
            initial.y = y0; initial.u = u0;
            [ initial , zeta0 ] = obj.mpc.get_zeta( initial );
            
            % initialize results struct
            results = struct;
            results.T = [ 0 ];
            results.U = [ u0( end , : ) ];
            results.Y = [ y0( end , : ) ];
            results.K = [ 0 ];
            results.R = [ ref(1,:) ];
            results.X = [ x0( end , : ) ];
            results.Z = []; % lifted states
            
            k = 1;
            while k < size( ref , 1 )
                
                % current time
                t = k * obj.mpc.params.Ts;
                
                % get current state and input with delays
                if k == 1
                    current.y = y0 * obj.scaledown.y;
                    current.u = u0 * obj.scaledown.u;
                elseif k < nd + 1
                    y = [ y0( k : end-1 , : ) ; results.Y ];
                    u = [ u0( k : end-1 , : ) ; results.U ];
                    current.y = y * obj.scaledown.y;
                    current.u = u * obj.scaledown.u;
                else
                    y = results.Y( end - nd : end , : );
                    u = results.U( end - nd : end , : );
                    current.y = y * obj.scaledown.y;
                    current.u = u * obj.scaledown.u;
                end
                
                % isolate the reference trajectory over the horizon
                if k + Np <= size( ref_sc , 1 )
                    refhor = ref_sc( k : k + Np , :);
                else
                    refhor = ref_sc( k : end , : );
                end
                    
%                 % isolate the shape_bounds over the horizon (TODO)
%                 if k + Np <= size( shape_bounds , 1 )
%                     shapehor = shape_bounds( k : k + Np , :);
%                 else
%                     shapehor = shape_bounds( k : end , : );
%                 end

                % get optimal input over horizon
                [ U , z ] = obj.mpc.get_mpcInput( current , refhor , shape_bounds );
                
                % isolate input for this step
                u_k_sc = U( 2 , : )';
                
                % scaleup the input for the system simulation
                u_k = obj.scaleup.u * u_k_sc;
                
                % simulate the system over one time-step
                x_k = results.X( end , : )';
                x_kp1 = obj.sys.simulate_Ts( x_k , u_k );
                y_kp1 = obj.sys.get_y( x_kp1' );
                
                % record updated results
                results.T = [ results.T ; t ];
                results.U = [ results.U ; u_k' ];
                results.Y = [ results.Y ; y_kp1 ];
                results.K = [ results.K ; k ];
                results.R = [ results.R ; ref( k , : ) ];
                results.X = [ results.X ; x_kp1' ];
                results.Z = [ results.Z ; z' ]; % current lifted state
                
                k = k + 1;  % increment step counter
            end
        end
    end
end

