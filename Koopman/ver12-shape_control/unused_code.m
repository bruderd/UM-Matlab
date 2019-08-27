% Old or unused versions of functions

%% sysid class

% get_scale (scale sim/exp data to be in range [-1 , 1])
function [ data_scaled , obj ] = get_scale( obj , data )
%scale: Scale sim/exp data to be in range [-1 , 1]
%    Also creates scaleup/scaledown matrices and saves as params
%    data - struct containing fields t , y , u (at least)
%    data_scaled - struct containing t , y , u , x (optional)

% get mean value in each dimension
y_mean = mean( data.y );
u_mean = mean( data.u );

% subtract the mean value to center data about zero
y_cent = data.y - y_mean;
u_cent = data.u - u_mean;

% get max absolute values in each dimension
y_maxabs = max( abs( y_cent ) );
u_maxabs = max( abs( u_cent ) );

% scale the data
data_scaled = struct;    % initialize
data_scaled.t = data.t;  % time is not scaled
data_scaled.y = y_cent ./ y_maxabs;
data_scaled.u = u_cent ./ u_maxabs;

% save scaling functions
y = sym( 'y' , [ 1 , obj.params.n ] );
u = sym( 'u' , [ 1 , obj.params.m ] );
y_scaledown = ( y - y_mean ) ./ y_maxabs;
u_scaledown = ( u - u_mean ) ./ u_maxabs;
obj.scaledown.y = matlabFunction( y_scaledown , 'Vars' , {y} );
obj.scaledown.u = matlabFunction( u_scaledown , 'Vars' , {u} );

y_scaleup = ( y .* y_maxabs ) + y_mean;
u_scaleup = ( u .* u_maxabs ) + u_mean;
obj.scaleup.y = matlabFunction( y_scaleup , 'Vars' , {y} );
obj.scaleup.u = matlabFunction( u_scaleup , 'Vars' , {u} );

%EVEN OLDER VERSION
%             % save the scaling matrices (note, these are meant to premultiply column vectors or postmultiply row vectors)
%             obj.params.scaleup.y = diag( y_maxabs );
%             obj.params.scaleup.u = diag( u_maxabs );
%             obj.params.scaledown.y = diag( y_maxabs .^ (-1) );
%             obj.params.scaledown.u = diag( u_maxabs .^ (-1) );
%
%             % do same for x if it is part of data struct
%             if ismember( 'x' , fields(data) )
%                 x_maxabs = max( abs( data.x ) );
%                 data_scaled.x = data.x ./ x_maxabs;
%                 obj.params.scaleup.x = diag( x_maxabs );
%                 obj.params.scaledown.x = diag( x_maxabs .^ (-1) );
%             end
%
%             % define zeta scaling matrix
%             scaledown_ydelays = kron( eye(obj.delays+1) , obj.params.scaledown.y );
%             scaledown_udelays = kron( eye(obj.delays) , obj.params.scaledown.u );
%             obj.params.scaledown.zeta = blkdiag( scaledown_ydelays , scaledown_udelays );
%             scaleup_ydelays = kron( eye(obj.delays+1) , obj.params.scaleup.y );
%             scaleup_udelays = kron( eye(obj.delays) , obj.params.scaleup.u );
%             obj.params.scaleup.zeta = blkdiag( scaleup_ydelays , scaleup_udelays );
%
%             % define z (lifted state) scaling matrix (THIS DOES NOT WORK FOR BASIS FUNCTIONS WITH SINES ANS COSINES)
%             zeta_scaledown = obj.params.scaledown.zeta * ones( obj.params.nzeta , 1);  % scaledown identity zeta vector
%             zeta_scaleup = obj.params.scaleup.zeta * ones( obj.params.nzeta , 1);  % scaleup identity zeta vector
%             zscales_down = obj.lift.full( zeta_scaledown ); % vector of scaling down factors
%             zscales_up = obj.lift.full( zeta_scaleup );     % vector of scaling up factors
%             obj.params.scaledown.z = diag( zscales_down );  % turn vector into scaling mtx
%             obj.params.scaleup.z = diag( zscales_up );     % turn vector into scaling mtx
end