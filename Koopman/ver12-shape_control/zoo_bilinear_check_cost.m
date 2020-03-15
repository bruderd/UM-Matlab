% zoo_bilinear_check_cost
%
% Moves along an actual robot trajectory, at each point, it calculates the
% MPC pseudoinput/input that should recreate the actual trajectory over the
% horizon, and saves the value of the cost function. It also evaluates the
% cost function for the actual inputs used by the system.
%
% This can be used to verify that the choice of MPC inputs is optimal with
% respect to the given model, it just may not have enough control authority
% to actually track the desired output.


function [ real , sim ] = zoo_bilinear_check_cost( trial , mpc )
%   trial - simulated trial with fields t, y, u
%   mpc - is an mpc class object

Nsteps = size( trial.t , 1 );

% outputs
real.cost = zeros( Nsteps , 1 ); % init
real.u = mpc.scaledown.u( trial.u );
real.nu = zeros( Nsteps , mpc.params.m );    % init
real.nu(1,:) = real.u(1,:); % same at first time step
real.t = trial.t;
real.y = zeros( Nsteps , mpc.params.n );

sim.cost = zeros( Nsteps , 1 ); % init
sim.nu = zeros( Nsteps , mpc.params.m );    % init
sim.u = zeros( Nsteps , mpc.params.m );     % init
sim.u(1,:) = real.u(1,:);
sim.nu(1,:) = sim.u(1,:);
sim.t = trial.t;

% simulate using the discrete Koopman model and the inputs from trial
zx = mpc.lift.zx( mpc.scaledown.y( trial.y( 1 , : ) )' );
real.zx(1,:) = zx';
real.y(1,:) = trial.y( 1 , : );
for i = 1 : Nsteps - 1
    zx_stack = kron( eye(length(mpc.params.u)) , zx );
    zxp1 = mpc.model.A * zx + mpc.model.B * zx_stack * real.u(i,:)';
    real.y(i+1,:) =  ( mpc.model.C * zxp1 )';
    real.zx(i+1,:) = zxp1;
    zx = zxp1;
end
sim.y = real.y; % init


Np = mpc.horizon;    % length of horizon
% step through trajectory and do calculations
for i = 1 : Nsteps - 1
    
    % get current state and input (no delays allowed)
    current.y = real.y( i , : );
    current.u = real.u( i , : );
%     current.zx = mpc.lift.zx( current.y' );
    current.zx = real.zx( i , : )';

    % isolate the reference trajectory over the horizon
    if i + Np <= size( real.t , 1 )
%         refhor = real.y( i : i + Np , :);
        refhor = real.zx( i : i + Np , :);
    else
%         refhor = real.y( i : end , : );     % repeat last entry
        refhor = real.zx( i : end , : );     % repeat last entry
    end

    % lift the reference trajectory (try to vectorize by modifying the zx lifting function
    ref_zx = zeros( size(refhor,1) , length(current.zx) );
    for j = 1 : size( refhor , 1 )
%         ref_zx(j,:) = mpc.lift.zx( refhor(j,:)' )';
        ref_zx(j,:) = refhor(j,:);
    end
    
    % vectorize the reference trajectory
    Yr = reshape( ref_zx' , [ ( Np + 1 ) * size(ref_zx,2) , 1 ] );
    
    % get optimal pseudoinput over horizon
    [ E , Zx ] = mpc.get_mpcInput_bilinear( current , refhor , [] );
    sim.nu(i+1,:) = E(2,:);     % isolate the pseudoinput for the next step
    Evec = reshape( E' , [ Np * mpc.params.m , 1 ] ); % stack into a vector
    
    % convert mpc pseudoinput to actual input
    zx_stack = kron( eye(length(mpc.params.u)) , Zx(2,:)' );
    zx0_stack = kron( eye(length(mpc.params.u)) , Zx(1,:)' );
    Bzx = mpc.model.B * zx_stack;
    Bzx0 = mpc.model.B * zx0_stack;
    sim.u(i+1,:) = lsqminnorm( Bzx , Bzx0 * sim.nu(i+1,:)' );
    
    % get cost matrices used by gurobi solver
    H = mpc.cost.get_H( current.y' );      
    G = mpc.cost.get_G( current.y' );
    D = mpc.cost.get_D( current.y' );
    f = ( current.zx' * G + Yr' * D )';
    
    % calculate cost of the mpc pseudoinputs vs the real pseudoinputs
    sim.cost(i+1,:) = Evec' * H * Evec + f' * Evec;
    
    
    % convert actual input to pseudoinput over horizon
    real_nuhor = zeros( Np , mpc.params.m );
    for k = 1 : Np
        zx0 = current.zx;
%         zxnow = mpc.lift.zx( real.y(i-1+k,:)' );
        zxnow = real.zx( i-1+k , : )';
        zx_stack = kron( eye(length(mpc.params.u)) , zxnow );
        zx0_stack = kron( eye(length(mpc.params.u)) , zx0 );
        Bzx = mpc.model.B * zx_stack;
        Bzx0 = mpc.model.B * zx0_stack;
        real_nuhor(k,:) = lsqminnorm( Bzx0 , Bzx * real.u(i-1+k,:)' );
    end
    real.nu(i+1,:) = real_nuhor(2,:);
    NUvec = reshape( real_nuhor' , [ Np * mpc.params.m , 1 ] ); % stack into a vector
    real.cost(i+1,:) = NUvec' * H * NUvec + f' * NUvec;
    
end


















