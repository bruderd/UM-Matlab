function unext = get_MPCinput_fb( know , xnow , xlast , ulast , load , ref , mpc , model)
% function to determine input at each time step. 
%   the "fb" at the end indicates the new feedback input "ulast"
%   NOTE: only compatible with models that have 0-1 delays.
%   ulast should be p x 2. First column is u_{k-2}, second is u_{k-1}

Np = mpc.params.Np; % number of steps in horizon
nd = model.params.nd;   % number of delays in model

% construct zeta
if nd == 0
    zeta = xnow;
else
    zeta = [ xnow ;...
             xlast ;...
             ulast(:,2) ];
end
    
% lift the state (must first identify lifting function corresponding to current model)
cd( ['..' , filesep , 'liftingFunctions'] );
lift = str2func( [ 'lift_' , model.params.systemName ] );
cd(['..' , filesep , 'mpc']);
z = lift( zeta, load );   % lift zeta and the estimate of the load

% isolate ref over horizon, if horizon goes beyond ref set extra points to the last point in ref
if (know + Np) <= size(ref,1) % this version is current step plus horizon
    Yr = reshape( ref( know : know+Np , : )' , [ model.params.ny * (Np+1) , 1 ] );
elseif (know + 1) <= size(ref,1)
    Yr = kron( ones(Np+1,1) , ref(end,:)' );  % vector where the last point in ref is repeated over entire horizon
    Yr(1 : model.params.ny * size(ref( know : end , : ) , 1)) = reshape( ref( know : end , : )' , [ model.params.ny * size( ref( know : end , : ) , 1) , 1 ] );    % replace first elements with remainder of ref trajectory
else
    Yr = kron( ones(Np+1,1) , ref(end,:)' );  % vector where the last point in ref is repeated over entire horizon
end


H = mpc.H;      % removed factor of 2 on 12/10/2018
f = ( z' * mpc.G + Yr' * mpc.D )';
A = mpc.L;
b = - mpc.M * z + mpc.c;

% Tack on "memory" constraint for last 2 inputs! (see page 183 in lab notebook)
Atack_top = [ [ speye( model.params.p ) ; -speye( model.params.p ) ] , sparse( 2*model.params.p , size(A,2) - model.params.p ) ];
Atack_bot = [ sparse( 2*model.params.p , model.params.p) , [ speye( model.params.p ) ; -speye( model.params.p ) ] , sparse( 2*model.params.p , size(A,2) - 2*model.params.p ) ];
Atack = [ Atack_top ; Atack_bot ];
btack = [ ulast(:,1) ; -ulast(:,1) ; ulast(:,2) ; -ulast(:,2) ]; %+ 1e-4; % added tolerance of 1e-4
A = [A ; Atack];    % tack on memory constraint
b = [b ; btack];

U = quadprog_gurobi( H , f , A , b );   % solve using gurobi (returns NaNs of cannot be solved)
% U = quadprog( H , f , A , b );     % solve using matlab

% set output
if all( isnan( U ) )    % if the problem was infeasable, hold same pressure as last time
    unext = ulast;
else
    unext = U( 2*model.params.p + 1 : 3 * model.params.p );    % the actual output of this function (third element of U)
end

end