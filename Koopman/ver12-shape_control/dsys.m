classdef dsys
    %dsys: Class with system properties and equations of motion inside
    %   This is for a generic system. Just change vf
    
    properties
        params struct;  % model parameters
    end
    
    methods
        function obj = dsys( params , varargin )
            %CLASS CONSTRUCTOR
            
            % set default parameter values
            params.sysName = 'bilinear-motor';
            params.n = 2;       % state dimension
            params.m = 1;       % input dimension
            params.Ts = 0.01;   % sampling time
            params.umax = 1;    % maxabs input value
            % other needed parameters (try to remove extranneous ones)
            params.nx = params.n;   % dimension of the full state (joint angles and joing velocities)
            params.ny = 2;   % dimension of measured output (mocap coordinates + end effector orientation)
            params.nu = params.m;  % dimension of the input (reference angle at each joint)
            obj.params = params;
            
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
        
        % saveclass
        function saveclass(obj)
            % save this system for later use
            dirname = [ 'systems' , filesep , obj.params.sysName ];
            unique_dirname = auto_rename( dirname , '(0)' );
            obj.params.sysName = erase( unique_dirname , ['systems', filesep] ) ;    % modify the system name parameter
            
            % create directory for the system, and save the class
            mkdir( unique_dirname );
            mkdir( [ unique_dirname , filesep , 'simulations' ] );  % make simulation subfolder
            fname = [ unique_dirname , filesep , obj.params.sysName, '.mat' ];
            save( fname , 'obj' );
        end
        
        %% dynamics
        
        % vf (dynamics as vector field)
        function f = vf( obj , t , x , u )
            %vf: Explicit dynamics for system
            
%             % Bilinear motor from Korda/Mezic MPC paper
%             La = 0.314;
%             Ra = 12.345;
%             km = 0.253;
%             J = 0.00441;
%             B = 0.00732;
%             tau1 = 1.47;
%             ua = 60;
%             
%             f = zeros(2,1);
%             f(1) = -( Ra / La ) * x(1) - ( km / La ) * x(2) * u(1) + ua / La;
%             f(2) = -( B / J ) * x(2) + ( km / J ) * x(1) * u(1) - tau1 / J; 
            
            % simple bilinear system
            f = zeros(2,1);
            f(1) = x(2);
            f(2) = -0.5*x(2) + x(1) * u(1) + 0.1 * u(1);
        end
        
        %% sensing
        
        % get_y (extracts the measured output from the full state)
        function y = get_y( obj , x )
            %get_y: Gets the output (in this case marker positions
            % vectorized) from the state (in this case Alpha)
            %   x - array containing one or more state values. Each row
            %    should be a state.
            %   y - array containing the corresponding output values. Each
            %    row is an output.
            
            y = x;
%             y = x(:,2);
        end
            
        %% simulation
        
        % simulate system under random "ramp and hold" inputs
        function sim = simulate( obj , tf , Tramp , varargin)
            %simulate_arm: simulate system under random "ramp and hold" inputs
            %   tf - length of simulation(s)
            %   Tramp - ramp period length
            %   Alpha - joint angles and velocities at each timestep
            %   markers - marker position at each time step [x0,y0,...,xn,yn ; ...]
            %   varargin - save on? (true/false)
            
            % save simulation results? Default is no
            if length(varargin) == 1
                saveon = varargin{1};
            else
                saveon = false;
            end
            
            % time steps
            tsteps = ( 0 : obj.params.Ts : tf )';    % all timesteps
            tswitch = ( 0 : Tramp : tf )';  % input switching points
            
            % table of inputs
            numPeriods = ceil( length(tswitch) / 2 );
            inputs_nohold = obj.params.umax .* ( 2*rand( numPeriods , obj.params.m ) - 1 );  % table of random inputs
            inputs_hold = reshape([inputs_nohold(:) inputs_nohold(:)]',2*size(inputs_nohold,1), []); % repeats rows of inputs so that we get a hold between ramps
            u = interp1( tswitch , inputs_hold( 1:length(tswitch) , : ) , tsteps );
            
            % initial condition (resting)
            x0 = zeros( obj.params.n , 1 );
            
            % simulate system
            [ t , x ] = ode45( @(t,x) obj.vf( t , x , u(floor(t/obj.params.Ts) + 1 , : )' ) , tsteps , x0);  %variable time step
%             x = ode5( @(t,x) obj.vf( t , x , u( floor(t/obj.params.Ts) + 1 , : )' ) , tsteps , x0);    % fixed time step 
%             t = tsteps;
            
            % define output
            sim.t = t;  % time vector
            sim.x = x;  % internal state of the system
            sim.y = obj.get_y( x );    % measured output
            sim.u = u;  % input
            sim.params = obj.params;    % parameters associated with the system
            
            % save results
            if saveon
                fname = [ 'systems' , filesep , obj.params.sysName , filesep , 'simulations' , filesep , 'tf-', num2str(tf) , 's_ramp-' , num2str(Tramp) , 's.mat' ];
                unique_fname = auto_rename( fname , '(0)' );
                save( unique_fname , '-struct' ,'sim' );
            end
        end
           
        % simulate_Ts (simulate system over a single time step)
        function [ x_kp1 ] = simulate_Ts( obj , x_k , u_k )
            %simulate_Ts: Simulate system over a single time step
            %   x_k - current value of state (full state Alpha = [alpha ; alphadot])
            %   u_k - input over next time step
            
            tstep = [ 0 , obj.params.Ts ];
            
            % simulate system
            [ ~ , x ] = ode45( @(t,x) obj.vf( t , x , u_k ) , tstep , x_k);    % with mass matrix, variable time step
            
            % set output, the state after one time step
            x_kp1 = x(end,:)';
        end
        
    end
end