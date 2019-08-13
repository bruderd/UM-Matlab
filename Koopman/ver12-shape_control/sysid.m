classdef sysid
    %sysid: Class for doing Koopman-based sysid
    %   Constructed from a sample set of data from an experiment or simulation.
    
    properties
        params struct;  % sysid parameters
        lift struct;    % contains matlab generated lifting functions
        basis struct;   % contains the observables for the system
        model struct;   % contains lifted system model of system
    end
    
    methods
        function obj = sysid( data , varargin )
            %sysid: Construct an instance of this class
            %   data - struct with fields t,x,u,y. Contains results of a
            %   simuation or experiment of the system to be identified
            %   (note: this is not necessarily the same data to be used for
            %   sysid, it just provides the dimensions of the state, input,
            %   etc.)
            
            obj.params = struct;    % initialize params struct
            obj.params.n = size( data.y , 2 );  % dimension of measured state
            obj.params.m = size( data.u , 2 );  % dimension of input
            obj.params.Ts = mean( data.t(2:end) - data.t(1:end-1) );    % sampling time
            obj.params.nd = 1;  % number of delays (defaults to 1) 
            
            % load in system parameters if any are provided
            if length(varargin) == 1
                obj.params.sysParams = varargin{1};
                obj.params.isfake = 1; % indicates whether class is built from a simulated system or real data
            else
                obj.params.isfake = 0;
            end
            
            obj.lift = struct;  % initialize lift struct
            obj.basis = struct; % initialize basis struct
            obj.model = struct; % initialize model struct
        end
        
        %% operations on simulation/experimental data
        
        % resample (resamples data with a desired time step)
        function data_resampled = resample( obj , data , Ts )
            % get query points
            tq = ( data.t(1) : Ts : data.t(end) )';
            
            data_resampled.t = tq;
            data_resampled.x = interp1( data.t , data.x , tq );
            data_resampled.u = interp1( data.t , data.u , tq );
            data_resampled.y = interp1( data.t , data.y , tq );
        end
        
        %% save the class object
        
        % save_class
        function obj = save_class( obj , varargin)
            %save_class: Saves the class as a .mat file
            %   If class is from a simulated system, it is saved in the
            %   subfolder corresponding to that system.
            %   If class if from a real system, it is saved in the generic
            %   folder /systems/fromData/.
            %   varargin = isupdate - indicates whether this is an update of a
            %   previously saved class (1) or a new class (0).
            
            % assume this is a new class file (not update) by default
            if length(varargin) == 1
                isupdate = varargin{1};
            else
                isupdate = 0;
            end
            
            % set the class name based on its parameters
            if isupdate
                classname = obj.params.classname;
            else
                dateString = datestr(now , 'yyyy-mm-dd_HH-MM'); % current date/time
                classname = [ 'n-' , num2str( obj.params.n ) , '_m-' , num2str( obj.params.m ) , '_del-' , num2str( obj.params.nd ) , '_' , dateString ];
                obj.params.classname = classname;   % create classname parameter
            end
    
            % save the class file
            model = obj;
            if obj.params.isfake    % check if class was created from simulated system
                dirname = [ 'systems' , filesep , obj.params.sysParams.sysName , filesep , 'models'];
                if ~isupdate
                    mkdir( dirname );
                end
                fname = [ dirname , filesep , classname, '.mat' ];
                save( fname , 'model' );
            else
                dirname = [ 'systems' , filesep , 'fromData' ];
                fname = [ dirname , filesep , classname, '.mat' ];
                save( fname , 'model' );
            end
        end
        
        %% defining observables
        
        % def_observables (define the basis of observable functions)
        function obj = def_observables( obj , type , degree )
            % def_observables: Defines the set of nonlinear observable
            % functions that will act as basis of Koopman subspace
            %   type - Type of functions, [cell array of strings].
            %       'armshape' - coefficients of shape polynomial
            %       'poly'  - monomials
            %       ... more to be added over time
            %   degree - Maximum degree/complexity of functions, [vector].
            
            % check that the type and degree have same dimension
            if size(type) ~= size(degree)
                error('inputs must be of the same size');
            end
            
            % define the low dimensional state with delays, called zeta
            x = sym('x', [obj.params.n, 1] , 'real');   % state variable x
            xd = sym('xd', [obj.params.nd * obj.params.n, 1] , 'real');   % state delays i.e. for 2 delays: [x_i-1, x_i-2]'
            ud = sym('ud', [obj.params.nd * obj.params.m, 1] , 'real');   % input delays i.e. for 2 delays: [u_i-1, u_i-2]'
            zeta = [x ; xd; ud];    % state variable with delays
            u = sym('u', [obj.params.m, 1] , 'real');   % input vector u
            obj.params.zeta = zeta;
            
            % construct the observables
            fullBasis = [];
            for i = 1 : length(degree)
                if strcmp( type{i} , 'armshape' )
                    obj = obj.def_armshapeLift( degree(i) );
                    fullBasis = [ fullBasis ; obj.basis.armshape ];
                elseif strcmp( type{i} , 'poly' )
                    obj = obj.def_polyLift( degree(i) );
                    fullBasis = [ fullBasis ; obj.basis.poly ];
                end
            end
            
            % define outputs
            obj.basis.full = fullBasis;
            obj.lift.full = matlabFunction( fullBasis , 'Vars' , {zeta} );
            obj.params.N = length( fullBasis ); % the dimension of the lefted state
            
%             % save the sysid class object (may not want to include this here, better done externally)
%             obj = obj.save_class( 0 );
        end
        
        % def_polyLift (defines polynomial basis functions)
        function obj = def_polyLift( obj , degree )
            %def_polyLift: Defines the lifting function that lifts state variable x to
            % space spanned by monomials with total degree less than or equal to
            % max_degree.
            %   e.g. 1 x1 x2 x1^2 x1x2 x2^2 ...
            
            zeta = obj.params.zeta; % get the symbolic unlifted state
            nzeta = length(zeta);
            maxDegree = degree;
            
            % Number of mononials, i.e. dimenstion of p(x)
            N = factorial(nzeta + maxDegree) / ( factorial(nzeta) * factorial(maxDegree) );
            
            % matrix of exponents (N x naug). Each row gives exponents for 1 monomial
            exponents = [];
            for i = 1:maxDegree
                exponents = [exponents; partitions(i, ones(1,nzeta))];
            end
            exponents = [exponents ; zeros(1,nzeta)];   % put constant at end of basis so state can be the first nzeta elements
            
            % create vector of orderd monomials (column vector)
            for i = 1:N
                polyBasis(i,1) = obj.get_monomial(zeta, exponents(i,:));
            end
            
%             % define matrix of exponents: columns=monomial term, rows=dimension of x
%             psi = exponents';
            
            % create the lifting function: zeta -> p(zeta)
            obj.lift.poly = matlabFunction(polyBasis, 'Vars', {zeta});
%             
%             % define derivative of lifted state with respect to x
%             dlift = jacobian(polyBasis,x);
%             matlabFunction(dlift, 'File', 'jacobianLift', 'Vars', {zeta});
            
            % output variables
            obj.basis.poly = polyBasis;    % symbolic vector of basis monomials, p(x)
%             ams.jacobianBasis = dlift;
%             params.N = N;   % dimension of polyBasis (including the state itself)
%             params.Np = N + p; % dimension of the lifted state
%             params.psi = psi;   % monomial exponent index function
%             params.x = x;   % symbolic state variable
%             params.u = u;   % symbolic input variable
%             params.xd = xd; % symbolic state delays
%             params.ud = ud; % symbolic input delays
            
        end
        
        % get_monomial (builds a monomial)
        function [ monomial ] = get_monomial( obj, x, exponents )
            %get_monomial: builds a monomial from symbolic vector x and a vector of
            %exponents
            %   e.g. x = [x1 x2]; exponents = [1 2]; =>  monomial = x1^1 * x2^2
            
            n = length(x);
            
            monomial = x(1)^exponents(1);
            for i = 2:n
                monomial = monomial * x(i)^exponents(i);
            end
        end
            
        % def_armshapeLift (defines the basis functions that describe the shape of robot arm)
        function obj = def_armshapeLift( obj , degree )
            
            % 2D case (will deal with the 3D case later)
            if ( ( obj.params.n-2 ) / obj.params.sysParams.Nmods ) == 2
                
                % define symbolic parameters
                zeta = obj.params.zeta;
                x = zeta( 1 : 2 : obj.params.n - 3 );   % x-coordinates of markers
                y = zeta( 2 : 2 : obj.params.n - 2 );    % y-coordinates of markers
                orient = zeta( obj.params.n-1 : obj.params.n );   % unit vector pointin in direction of end effector
                
                % define intermediate variables
                points = [ x , y ]; % matrix of marker coordinates. Each row is coords of a marker
                positions = obj.params.sysParams.markerPos(2:end); % row vector of marker positions [0,1]
                
                % generate virtual points to provide slope constraint at the base & end
                startpoint = [ 0 , 1e-2 ];
                endpoint = ( orient' )*1e-2 + points(end,:);    % add point extended from last link
                points_supp = [0 , 0 ; startpoint ; points ; endpoint];
                
                % generate A matrix for least squares problem
                positions_supp = [ 0 , 1e-2 , positions , 1+1e-2 ];
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
                
                % define output
                armshapeBasis = [coeffs_vec_x ; coeffs_vec_y ];
                obj.basis.armshape = armshapeBasis;
                
                % create the lifting function: zeta -> p(zeta)
                obj.lift.armshape = matlabFunction(armshapeBasis, 'Vars', {zeta});
            else
                error('Number of marker postions must be number of modules (for now)');
            end
        end
        
        %% fitting Koopman operator and A,B,C system matrices
        
        % get_zeta (adds a zeta field to a test data struct)
        function [ data_out , zeta ] = get_zeta( obj , data_in )
            %get_zeta: Adds a zeta field to a test data struct
            %   data_in - struct with t , x , y , u fields
            %   zeta - [ y , yd1 , yd2 , ... , ud1 , ud2 , ... ]
            
            data_out = data_in;
            
            % add the zeta field
            for i = obj.params.nd + 1 : length( data_in.t )
                ind = i - obj.params.nd;    % current timestep index
                y = data_in.y( i , : );
                ydel = zeros( 1 , obj.params.nd * obj.params.n );
                udel = zeros( 1 , obj.params.nd * obj.params.m );
                for j = 1 : obj.params.nd
                    fillrange_y = obj.params.n * (j - 1) + 1 : obj.params.n * j;
                    fillrange_u = obj.params.m * (j - 1) + 1 : obj.params.m * j;
                    ydel(1 , fillrange_y) = data_in.y( i - j , : );
                    udel(1 , fillrange_u) = data_in.u( i - j , : );
                end
                zetak = [ y , ydel , udel ];
                data_out.zeta( ind , : ) = zetak;
                data_out.uzeta( ind , : ) = data_in.u( i , : );    % current timestep with zeta (input starting at current timestep)
            end
            zeta = data_out.zeta;
        end
        
        % get_snapshotPairs (convert time-series data into snapshot pairs)
        function snapshotPairs = get_snapshotPairs( obj , data , varargin )
            %get_snapshotPairs: Convert time-series data into a set of num
            %snapshot pairs.
            %   data - struct with fields x , y , u , t , zeta
            %   varargin = num - number of snapshot pairs to be taken
            
            % set the number of snapshot pairs to be taken
            num_max = size( data.zeta , 1 ) - 1; % maximum number of snapshot pairs
            if length(varargin) == 1
                num = varargin{1};
                if num > num_max - 1
                    message = [ 'Number of snapshot pairs cannot exceed ' , num2str(num_max) , '. Taking ' , num2str(num_max) , ' pairs instead.' ];
                    disp(message);
                    num = num_max;
                end
            else
                num = num_max;
            end
            
            % separate data into 'before' and 'after' time step
            before.zeta = data.zeta( 1:end-1 , : );
            after.zeta = data.zeta( 2:end , : );
            u = data.uzeta( 1:end-1 , : );    % input that happens between before.zeta and after.zeta
            
            % randomly select num snapshot pairs
            total = num_max;
            s = RandStream('mlfg6331_64'); 
            index = datasample(s , 1:total, num , 'Replace' , false);
            
            snapshotPairs.alpha = before.zeta( index , : ); 
            snapshotPairs.beta = after.zeta( index , : );
            snapshotPairs.u = u( index , : );
        end
        
        % get_Koopman (Find the best possible koopman operator from snapshot pairs)
        function [ U , koopData ] = get_Koopman( obj ,  snapshotPairs , varargin )
            %get_KoopmanConstGen: Find the best possible koopman operator given
            %snapshot pairs using constraint generation to deal with large data sets.
            %   varargin = lasso weighting parameter. lasso >> 1 approaches least squares solution 
            
            if length(varargin) == 1
                lasso = varargin{1};
            else
                lasso = 100 * obj.params.N;   % defualt value of the lasso parameter (should emulate least squares)
            end
            
            disp('Finding Koopman operator approximation...');
            
            % Extract snapshot pairs
            [x,y,u] = deal( snapshotPairs.alpha , snapshotPairs.beta , snapshotPairs.u );
            
            % Build matrices
            [n, m] = deal( obj.params.n , obj.params.m );
            Nm = obj.params.N + m;   % dimension of zeta plus input
            
            Px = zeros(length(x), Nm);
            Py = zeros(length(x), Nm);
            for i = 1:length(x)
                psix = obj.lift.full( x(i,:)' )';   
                psiy = obj.lift.full( y(i,:)' )';
                Px(i,:) = [ psix , u(i,:) ];
                Py(i,:) = [ psiy , zeros(1,m) ];     % exclude u from Py (could also use same u as Px)
            end
            
            % Call function that solves QP problem
            Uvec = obj.solve_KoopmanQP( Px , Py , lasso);
            U = reshape(Uvec, [Nm,Nm]); % Koopman operator matrix
            
            % other usefule outputs
            koopData.Px = Px( : , 1 : obj.params.N );   % only want state, not input
            koopData.Py = Py( : , 1 : obj.params.N );
            koopData.u = u;
            koopData.alpha = snapshotPairs.alpha;
        end
        
        % solve_KoopmanQP (Finds the Koopman operator using Lasso regression)
        function Uvec = solve_KoopmanQP( obj , Px , Py , lasso )
            %solve_KoopmanQP: Finds the Koopman operator for a given set of
            %data using the lasso regression method.
            %   Lasso method
            %   x is vectorized Koopman operator, decomposed into positive and negative parts of each entry x = [u11+, ..., uNN+, u11-, ... , uNN-]';
            %   Uvec = M * x, where M subtracts the + and - parts of each entry: uij+ - uij-

            nx = obj.params.N^2;
            Nm = obj.params.N + obj.params.m;
            
            M = [speye(Nm^2) , -speye(Nm^2)];
            
            PxTPx = Px' * Px;
            PxTPy = Px' * Py;
            ATA = kron(speye(Nm) , PxTPx);  % repeat blocks diagonally N times
            ATb = reshape(PxTPy, [Nm^2 , 1]);
            
            % L2 error as cost function
            preH = ATA * M;
            H = M' * preH;
            f = -M' * ATb;
            
            % L1 regularization enforced as constraint
            t = lasso;
            Aq = [ -speye(2*Nm^2) ; ones(1 , 2*Nm^2) ];
            bq = [ zeros(2*Nm^2 , 1) ; t ];
            
            % Solve the quadratic program
            [x , results] = quadprog_gurobi( H , f , Aq , bq );       % use gurobi to solve
            % options = optimoptions('quadprog', 'Display', 'iter');
            % [ x, fval, exitflag ] = quadprog(H, f, Aq, bq, [], [], [], [], [],options);      % use matlab to solve
            
            % Recover Uvec from the optimization variable
            xout = M * x;
            
            % Set output
            Uvec = xout;
        end
        
        % get_ABC (Extracts the A,B,C matrices from Koopman matrix)
        function [ out , obj ] = get_ABC( obj , U , koopData )
            %get_ABC: Defines the A, B, and C matrices that describe the linear
            %lifted system z+ = Az + Bu, x = Cz.
            %   out - struct with fields A, B, C, sys, params, ...
            %   obj.model - property of struct which stores the model
            
            UT = U';    % transpose of koopman operator
            
            A = UT( 1 : obj.params.N , 1 : obj.params.N );
            B = UT( 1 : obj.params.N , obj.params.N+1 : end );
            
            % C selects the first ny entries of the lifted state (so output can be different than state)
            Cy = [eye(obj.params.n), zeros(obj.params.n , obj.params.N - obj.params.n)];   

            % matrix to recover the state with delays, i.e. zeta
            nzeta = obj.params.n + obj.params.nd * ( obj.params.n + obj.params.m );
            Cz = [eye( nzeta ), zeros( nzeta , obj.params.N - nzeta)];
            
            % % solve for C matrix that best recovers state from lifted state
            % Ctranspose = koopData.Px \ koopData.x ;
            % C = Ctranspose' ;
            
            % find M matrix that (approximately) projects a lifted point onto the subset of all legitimate lifted points in R^N
            Px = koopData.Px; Py = koopData.Py;
            U = koopData.u;
            L = zeros( size(Px,1) , obj.params.N );
            for i = 1 : size( Px , 1)
                L(i,:) = ( A * Px(i,:)' + B * U(i,:)' )' ;        % with input
                %     L(i,:) = ( A * Px(i,:)' )' ;        % without input
            end
            R = zeros( size(L,1) , obj.params.N );
            for i = 1 : size( L , 1 )
                %     R(i,:) = ( stateLift( Cz * L(1,:)' ) )' ;
                R(i,:) = Py(i,:) ;  % debug. This should make M the identity
            end
            Mtranspose = L \ R;
            M = Mtranspose';
            
            % define outputs
            out.A = M*A;  % edited to include projection M, 12/11/2018
            out.B = M*B;  % edited to include projection M, 12/11/2018
            out.Asim = A;
            out.Bsim = B;
            out.C = Cy;
            out.Cz = Cz;
            out.M = M;
            out.sys = ss( out.A , out.B , Cy , 0 , obj.params.Ts );  % discrete state space system object
            out.params = obj.params;    % save system parameters as part of system struct    
            
            % add model substruct to the class
            obj.model = out;
        end
        
        %% validate performance of a fitted model
        
        % val_liftedSys (compares model simulation to real data)
        function results = val_liftedSys( obj , liftedSys , valdata )
            %val_liftedSys: Compares a model simulation to real data
            %   liftedSys - struct with fields A, B, C, sys, ...
            %   valdata - struct with fields t, y, u (at least)
            %   results - struct with simulation results and error calculations
            
            % shift real data so delays can be accounted for
            index0 = obj.params.nd + 1;  % index of the first state
            treal = valdata.t(index0 : end);    % start simulation late so delays can be taken into account
            yreal = valdata.y(index0 : end , :);
            ureal = valdata.u(index0 : end , :);
            
            % set initial condition
            valdata_wzeta = get_zeta( obj , valdata );
            zeta0 = valdata_wzeta.zeta(1,:)';    % initial state with delays
            z0 = obj.lift.full( zeta0 );    % initial lifted state
            
            % simulate lifted linear model
            [ ysim , tsim , zsim ] = lsim(liftedSys.sys, ureal , treal , z0);
            
            % save simulations in output struct
            results = struct;
            results.t = treal; 
            results.sim.t = tsim;
            results.sim.y = ysim;
            results.real.t = treal;
            results.real.y = yreal;
            
            % save error info
            results.error.abs = abs( ysim - yreal );
            results.error.mean = mean( results.error.abs , 1 );
            results.error.rmse = sqrt( sum( (ysim - yreal).^2 , 1 ) ./ length(treal) );
        end
        
        % plot_comparison (plots a comparison between simulation and real data)
        function plot_comparison( obj , simdata , realdata )
            %plot_comparison: plots a comparison between simulation and real data.
            
            fig = figure;   % create new figure
            for i = 1 : obj.params.n
                subplot( obj.params.n , 1 , i );
                hold on;
                plot( realdata.t , realdata.y( : , i ) , 'b' );
                plot( simdata.t , simdata.y( : , i ) , 'r' );
                hold off;
            end
        end
        
    end
end
























