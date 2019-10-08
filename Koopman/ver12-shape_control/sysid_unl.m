classdef sysid_unl < sysid
    %unl: sysid class with nonlinear input
    %   Detailed explanation goes here
    
    properties
        basis_type;     % types of nonlinear functions in e2u mapping
        basis_degree;   % degree of functions in e2u mapping
        
        e2u_basis;      % set of nonlinear functions in e2u mapping
        e2u_lift;       % function for evaluating basis functions
        
        e2u;            % struct to hold e2u function and nnetwork
    end
    
    methods
        function obj = sysid_unl( data4sysid , options , varargin )
            %CLASS CONSTRUCTOR
            %   Detailed explanation goes here
            
            % Create instance of sysid class
            obj = obj@sysid( data4sysid , varargin{:} );
            
            % Flag this as an object with a nonlinear input
            obj.params.NLinput = 1;
            
            % set default values of arguments that are specific to this class
            obj.basis_type = { 'poly' };
            obj.basis_degree = [ 1 ];
            
            % set values of arguments that are specific to this class
            obj = parse_args( obj , options{:} );
            
            % train model(s)
            obj = obj.train_models;
        end
        
        % parse_args: Parses the Name, Value pairs in varargin
        function obj = parse_args( obj , varargin )
            %parse_args: Parses the Name, Value pairs in varargin of the
            % constructor, and assigns property values
            for idx = 1:2:length(varargin)
                obj.(varargin{idx}) = varargin{idx+1} ;
            end
        end
        
        % get_e: calculate the error of the koopman model (w/o input) at each time step
        function [ e , data_out ] = get_e( obj , data_in )
            % get_e: calculate the error of the koopman model (w/o input) at each time step
            %   data_in - struct containing time series data or snapshot
            %     pairs
            
            if isfield( data_in , 't' )     % assume time series data
                if ~isfield( data_in , 'zeta' )
                    [ data_in , zeta ] = obj.get_zeta( data_in );
                end
                e = zeros( size( zeta , 1 ) - 1 , obj.params.N );
                for i = 1 : size( zeta , 1 ) - 1
                    e(i,:) = ( obj.lift.full( zeta(i+1,:)' ) - obj.model.A * obj.lift.full( zeta(i,:)' ) )';
                end
                data_out = data_in;
                data_out.e = e;
            else    % assume snapshot pairs
                e = zeros( size( data_in.alpha , 1 ) , obj.params.N );
                for i = 1 : size( data_in.alpha , 1 )
                    e(i,1) = ( data_in.beta(i,:)' - obj.model.A * data_in.alpha(i,:)' )';
                end
                data_out = data_in;
                data_out.e = e;
            end
        end
        
        % get_nu: calculate a reduced dimension version of e using pca
        function [ nu , Beta , data_out ] = get_nu( obj , data_in )
            
            if isfield( data_in , 'e' ) % assume data_in is struct containing time series or snapshot pairs 
                e = data_in.e;
            else    % assume data_in is just an e_matrix
                e = data_in;
            end
            
            % use pca to find most important components
            [ coeffs , ~ , ~ , ~ , explained , ~ ] = pca( e );
            
            % take enough components to explain > 90% of the data
            num_pcs = 1;
            while sum( explained(1:num_pcs) ) < 90
                num_pcs = num_pcs + 1;
            end
            
            % define projection matrix from e to nu
            Beta = coeffs( : , 1 : num_pcs )';
%             Beta = eye( size(e,2) );  % make nu same as e  % MUST UNDO THIS!!!!!!!!!!!!!!!!!!!!!!!
            
            % calculate nu for each value of e in data_in
            nu = e * Beta';
            
            % add nu as a field to data_in (if it is a struct)
            data_out = data_in; % just pass through
            if isfield( data_in , 'e' ) % assume data_in is struct containing time series or snapshot pairs 
                data_out.nu = nu;
            end
        end
        
        % def_e2u_basis (define the basis of observable functions) (UNUSED 2019-10-08)
        function obj = def_e2u_basis( obj , type , degree )
            % def_e2u_basis: Defines the set of nonlinear observable
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
            zeta_ = [x ; xd; ud];    % state variable with delays
            e = sym('e', [obj.params.N, 1] , 'real');   % linear input vector e
            zeta = [ e ; zeta_ ];   % this is argument into e2u. Only calling it zeta so it will work with other functions
            obj.params.zeta = zeta; % NOTE: this is overwriting the real zeta!!!
            
            % construct the observables
            fullBasis = zeta;  % first nzeta observables should always be the measured state and delays
            for i = 1 : length(degree)
                if strcmp( type{i} , 'armshape' )
                    [ ~ , armshape ] = obj.def_armshapeLift( degree(i) );
                    fullBasis = [ fullBasis ; armshape ];
                elseif strcmp( type{i} , 'poly' )
                    [ ~ , polyb ] = obj.def_polyLift( degree(i) );
                    fullBasis = [ fullBasis ; polyb( obj.params.N+obj.params.nzeta+1 : end ) ]; % don't include first nzeta states because they will be repeats
                elseif strcmp( type{i} , 'fourier' )
                    [ ~ , fourier ] = obj.def_fourierLift( degree(i) );
                    fullBasis = [ fullBasis ; fourier ];
                elseif strcmp( type{i} , 'fourier_sparser' )
                    [ ~ , fourier_sparser ] = obj.def_fourierLift_sparser( degree(i) );
                    fullBasis = [ fullBasis ; fourier_sparser ];
                elseif strcmp( type{i} , 'atan' )
                    [ ~ , atanb ] = obj.def_atanLift( degree(i) );
                    fullBasis = [ fullBasis ; atanb ];
                elseif strcmp( type{i} , 'gaussian' )
                    [ ~ , gaussian ] = obj.def_gaussianLift( degree(i) );
                    fullBasis = [ fullBasis ; gaussian ];
                elseif strcmp( type{i} , 'hermite' )
                    [ ~ , hermite ] = obj.def_hermiteLift( degree(i) );
                    fullBasis = [ fullBasis ; hermite ];
                end
            end
            
            % add a constant term to the end of the set
            fullBasis = [ fullBasis ; sym(1) ];
            
            % define outputs
            obj.e2u_basis = fullBasis;
            obj.e2u_lift = matlabFunction( fullBasis , 'Vars' , {zeta} );
        end
        
        
    end
end




















