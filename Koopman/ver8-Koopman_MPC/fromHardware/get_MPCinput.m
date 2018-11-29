function unext = get_MPCinput( know , xnow , ref , mpc , model)
% function to determine input at each time step

Np = mpc.params.Np; % number of steps in horizon

% need persistent variables in case there are delays
persistent k x u

% know = floor( tnow / mpc.Ts );    % current timestep

% initialize values of persistent variables
if isempty(k)
    k = know;
    x = kron( ones(1 , model.params.nd+1) , xnow );   % repeat initial condition for time steps that occur before k = 0
    u = kron( ones(1 , model.params.nd+1) , zeros(model.params.p , 1) );    % assume input is zero before k = 0
end

% update values of persistent variables
if know == k(end)       % still at same timestep?
    unext = u(: , end);     % u takes same value 
    return
else
    nd = model.params.nd;
    if nd == 0
        zeta = xnow;
    else
        zeta = [ xnow ;...
                reshape( x( : , end-nd+1 : end ) , [model.params.n * nd , 1] ) ;...
                reshape( u( : , end-nd+1 : end ) , [model.params.p * nd , 1] ) ];
    end
    
    % lift the state (must first identify lifting function corresponding to current model)
    cd( ['model' , filesep , 'liftFunctions'] );
    lift = str2func( [ 'lift_' , model.params.systemName ] );
    cd(['..' , filesep , '..']);
    z = lift(zeta);   % lift zeta
    
    % isolate ref over horizon, if horizon goes beyond ref set extra points to the last point in ref
    if (know + Np) <= size(ref,1) % this version is current step plus horizon
        Yr = reshape( ref( know : know+Np , : )' , [ model.params.ny * (Np+1) , 1 ] );
    elseif (know + 1) <= size(ref,1)
        Yr = kron( ones(Np+1,1) , ref(end,:)' );  % vector where the last point in ref is repeated over entire horizon
        Yr(1 : model.params.ny * size(ref( know : end , : ) , 1)) = reshape( ref( know : end , : )' , [ model.params.ny * size( ref( know : end , : ) , 1) , 1 ] );    % replace first elements with remainder of ref trajectory
    else
        Yr = kron( ones(Np+1,1) , ref(end,:)' );  % vector where the last point in ref is repeated over entire horizon
    end
%     if (know + Np) <= size(ref,1)     % this version is current step included in horizon
%         Yr = reshape( ref( know+1 : know+Np , : ) , [ model.params.ny * Np , 1 ] );
%     elseif (know + 1) <= size(ref,1)
%         Yr = kron( ones(Np,1) , ref(end,:)' );  % vector where the last point in ref is repeated over entire horizon
%         Yr(1 : model.params.ny * size(ref( know+1 : end , : ) , 1)) = reshape( ref( know+1 : end , : ) , [ model.params.ny * size( ref( know+1 : end , : ) , 1) , 1 ] );    % replace first elements with remainder of ref trajectory
%     else
%         Yr = kron( ones(Np,1) , ref(end,:)' );  % vector where the last point in ref is repeated over entire horizon
%     end
end

    H = 2 * mpc.H;
    f = ( z' * mpc.G + Yr' * mpc.D )';
    A = mpc.L;
    b = - mpc.M * z + mpc.c;
    U = quadprog_gurobi( H , f , A , b );   % solve using gurobi
%     U = quadprog( H , f , A , b );     % solve using matlab
    unext = U( 1 : model.params.p );    % the actual output of this function
    
    % store values of the state and input at the "beginning" of a timestep
    k = [k , know];     % vector of all previous timestamps
    x = [x , xnow];     % vector of all previous states
    u = [u , unext];
    
end


